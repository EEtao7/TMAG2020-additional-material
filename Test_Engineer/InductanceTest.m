classdef InductanceTest < TestProject
    properties
        % test items(测试项目)
        MyModel
        Inductance
        Ld
        Lq
        documentsaddress = 'E:\My_Study\ZJU_master_study\master research\Low speed servo motor\Matlab_and_Simulink_documents\figure documents\Inductance';
        figureaddress
        
        % adjust test parameters（可调节的测试参数）
        theta_min_deg
        theta_max_deg
        dtheta_deg
        
        Motor_flag
        Rotor_flag
        MagneticField_flag = 0
        Im = 3
        initial_angle
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
            Pi=pi;
            deg=Pi/180.;
            
            % some useful coefficients to transfer the values
            Periodic_coeff = 1; % when caculating the Inductance and Flux by femm Periodic BC, you need to multiplicate it. Not to the torque in femm.
            %             Park_coeff = 1.5; % when caculating the torque and power by equations in dq coordinate system.
            
            %% Other User-defined  Parameters
            % Set the mechanical angle alpha
            alpha_min = obj.theta_min_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            alpha_max = obj.theta_max_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            dalpha = obj.dtheta_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            
            %% Initialize the Variables
            obj.Inductance=zeros(obj.size_theta,1);
            
            %% Start a series of finite element caculations to Flux_PM of VPM
            if(obj.MagneticField_flag ==1)
                openfemm(0);
            elseif(obj.MagneticField_flag ==0)
                openfemm(1);
            end
            opendocument(obj.MyModel);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            mi_saveas('temp.fem');
            % Set the current in the circuit coils
            mi_setcurrent('A', obj.Im);
            mi_setcurrent('B', -0.5*obj.Im);
            mi_setcurrent('C', -0.5*obj.Im);
            % Set the Hc of PM material
            mi_modifymaterial('PM',3,0);
            
            for alpha_deg = alpha_min:dalpha:alpha_max
                m = round((alpha_deg-alpha_min)/dalpha)+1;
                % useful skills to caculate the time of circulation
                starttime = clock;
                
                %% Set the initial rotor position
                tta0 = obj.initial_angle/obj.testmotor.rotor.pole_of_pairs_Rotor; %make the d axis is coinside with A axis
                if (obj.Motor_flag==1)
                    ttad = -round(tta0*obj.testmotor.rotor.pole_of_pairs_Rotor);% rotor angle in degrees clockwise rotation
                elseif(obj.Motor_flag==2)
                    ttad = round(tta0*obj.testmotor.rotor.pole_of_pairs_Rotor); % rotor angle in degrees anticlockwise rotation
                end
                if (obj.Motor_flag==1)
                    ttaq = -round(90+obj.initial_angle); % rotor angle in degrees clockwise rotation
                elseif(obj.Motor_flag==2)
                    ttaq = round(90+obj.initial_angle); % rotor angle in degrees anticlockwise rotation
                end
                if (obj.Motor_flag==1)
                    Rotor_angle=-(m-1)*dalpha-tta0; % rotor angle in degrees VPM clockwise rotation
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
                circuit_B = mo_getcircuitproperties('B');
                obj.Inductance(m) = Periodic_coeff*0.6666667*((circuit_A(3)-circuit_B(3))/obj.Im);
                
                % save the Ld and Lq
                if (mod(round(Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor)-ttad, 360) == 0)
                    obj.Ld = obj.Inductance(m);
                elseif(mod(round(Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor)-ttaq, 360) == 0)
                    obj.Lq = obj.Inductance(m);
                end
                
                if(obj.MagneticField_flag ==1)
                    % save the magneticfield figure
                    mo_showdensityplot(0, 0, 3, 0, 'mag');
                    x1 = -105;
                    y1 = -105;
                    x2 = 105;
                    y2 = 105;
                    address = 'E:\My_Picture\';
                    name = 'test Motor';
                    number = m;
                    MagneticfieldFigure(x1, y1, x2, y2, address, name, number);
                end
                
                %% Close the Postprocessor and Delete the data of FEA
                mo_close;
                delete('temp.ans');
                fprintf('% i of % i :: %f seconds ::  %f H \n',round(Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor), obj.theta_max_deg, etime(clock,starttime), obj.Inductance(m));
                
            end
            %% Clean up after all finite element cacualation are finished
            closefemm;
            delete('temp.fem');
            fprintf('Ld = %f mH and Lq = %f mH ',obj.Ld, obj.Lq);
        end
        
        function Save_Data(obj)
            figure();
            plot((obj.theta_min_deg:obj.dtheta_deg:obj.theta_max_deg),obj.Inductance);
            xlabel('Electrical Angle, Degrees');
            ylabel('Inductance, H');
            title('Inductance in VPM');
            grid on;
            % save the figure document and useful parameters
            figurename = strcat(obj.testmotor.name,'_Inductance.fig');
            obj.figureaddress = strcat(obj.documentsaddress,'\', figurename);
            savefig(obj.figureaddress);
            obj.testmotor.controlparameters.Inductance =obj.Inductance;
            obj.testmotor.controlparameters.Ld =obj.Ld;
            obj.testmotor.controlparameters.Lq =obj.Lq;
        end
    end
end