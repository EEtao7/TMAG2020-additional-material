classdef FluxTest < TestProject
    properties
        MyModel
        flux_PM
        documentsaddress = 'E:\My_Study\ZJU_master_study\master research\Low speed servo motor\Matlab_and_Simulink_documents\figure documents\Flux linkage';
        figureaddress
        
        % adjust parameters
        Motor_flag
        Rotor_flag
        
        theta_min_deg
        theta_max_deg
        dtheta_deg
    end
    
    properties(Dependent)
        size_theta
    end
    methods
        function size_theta = get.size_theta(obj)
            size_theta = round((obj.theta_max_deg-obj.theta_min_deg)/obj.dtheta_deg)+1;
        end
        
        function Load_motor(obj)
            %% Load the Testmotor and Define some varieties
            obj.MyModel = obj.testmotor.fileaddress;
            
        end
        function FEA_analysis(obj)
            %% Useful Coeffiencys
            % helpful unit definitions
            %             Pi=pi;
            %             deg=Pi/180.;
            
            % some useful coefficients to transfer the values
            Periodic_coeff = 1; % when caculating the Inductance and Flux by femm Periodic BC, you need to multiplicate it. Not to the torque in femm.
            %             Park_coeff = 1.5; % when caculating the torque and power by equations in dq coordinate system.
            
            %% Other User-defined  Parameters
            % Set the mechanical angle alpha
            initial_angle = 0;
            alpha_min = obj.theta_min_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            alpha_max = obj.theta_max_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            dalpha = obj.dtheta_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            
            %% Initialize the Variables
            obj.flux_PM=zeros(obj.size_theta,1);
            
            %% Start a series of finite element caculations to Flux_PM of VPM
            openfemm(1);
            opendocument(obj.MyModel);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            mi_saveas('temp.fem');
            
            % Set the current in the circuit coils
            mi_setcurrent('A',0);
            mi_setcurrent('B',0);
            mi_setcurrent('C',0);
            
            for alpha_deg = alpha_min:dalpha:alpha_max
                m = round((alpha_deg-alpha_min)/dalpha)+1;
                % useful skills to caculate the time of circulation
                starttime = clock;
                
                %% Set the initial rotor position
                tta0 = initial_angle/obj.testmotor.rotor.pole_of_pairs_Rotor; %make the d axis is coinside with A axis
                if (obj.Motor_flag==1)
                    Rotor_angle=-(m-1)*dalpha+tta0; % rotor angle in degrees VPM clockwise rotation
                elseif(obj.Motor_flag==2)
                    Rotor_angle=(m-1)*dalpha+tta0; % rotor angle in degrees PMSM anticlockwise rotation
                end
                if(obj.Rotor_flag==1)
                    mi_modifyboundprop('AGE',11,Rotor_angle);% outer rotor
                elseif(obj.Rotor_flag==2)
                    mi_modifyboundprop('AGE',10,Rotor_angle);% inner rotor
                end
                
                %% Finite Elements Analysis Start
                mi_analyze(1);
                mi_loadsolution;
                mo_smooth('off');  % flux smoothing algorithm is off to increase the speed of code.
                
                %% Magnetic Postprocessor Command
                circuit_A = mo_getcircuitproperties('A');
                obj.flux_PM(m) = 1.5*Periodic_coeff*circuit_A(3);
                
                %% Close the Postprocessor and Delete the data of FEA
                mo_close;
                delete('temp.ans');
                fprintf('% i of % i :: %f seconds ::  %f Web \n',Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor,obj.theta_max_deg,etime(clock,starttime),obj.flux_PM(m));
                
            end
            %% Clean up after all finite element cacualation are finished
            closefemm;
            delete('temp.fem');
        end
        function Save_Data(obj)
%             plot((obj.theta_min_deg:obj.dtheta_deg:obj.theta_max_deg),obj.flux_PM);
%             xlabel('Electrical Angle, Degrees');
%             ylabel('Flux PM, Web');
%             title('Flux PM in VPM');
%             grid on;
%             % save the figure document and useful parameters
%             figurename = strcat(obj.testmotor.name,'_Flux_PM.fig');
%             obj.figureaddress = strcat(obj.documentsaddress,'\', figurename);
%             savefig(obj.figureaddress);
            obj.testmotor.controlparameters.flux_PM =obj.flux_PM;
        end
    end
end