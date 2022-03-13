classdef ReluctanceTorqueTest < TestProject
    properties
        % test items(测试项目)
        MyModel
        reluctance_torque
        documentsaddress = 'E:\My_Study\ZJU_master_study\master research\Maxium Efficiency Control of VPM\Low speed servo motor\Matlab_and_Simulink_documents\figure documents\Reluctance torque';
        figureaddress
        
        % adjust test parameters（可调节的测试参数）
        MyLowestHarmonic
        theta_min_deg
        theta_max_deg
        dtheta_deg
        
        Im = 1
        gamma_min_deg
        gamma_max_deg
        dgamma_deg
        
    end
    properties(Dependent)
        size_theta
        size_gamma
    end
    methods
        function size_theta = get.size_theta(obj)
            size_theta = round((obj.theta_max_deg-obj.theta_min_deg)/obj.dtheta_deg)+1;
        end
        function size_gamma = get.size_gamma(obj)
            size_gamma = round((obj.gamma_max_deg-obj.gamma_min_deg)/obj.dgamma_deg)+1;
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
            %             Periodic_coeff = 1; % when caculating the Inductance and Flux by femm Periodic BC, you need to multiplicate it. Not to the torque in femm.
            %             Park_coeff = 1.5; % when caculating the torque and power by equations in dq coordinate system.
            
            %% Other User-defined  Parameters
            % Set the mechanical angle alpha
            initial_angle = 0;
            obj.MyLowestHarmonic = 1;
            alpha_min = obj.theta_min_deg/obj.MyLowestHarmonic/obj.testmotor.rotor.pole_of_pairs_Rotor;
            alpha_max = obj.theta_max_deg/obj.MyLowestHarmonic/obj.testmotor.rotor.pole_of_pairs_Rotor;
            dalpha = obj.dtheta_deg/obj.MyLowestHarmonic/obj.testmotor.rotor.pole_of_pairs_Rotor;
            
            % Set the advanced current angle gamma
            gamma_min = obj.gamma_min_deg;
            gamma_max = obj.gamma_max_deg;
            dgamma = obj.dgamma_deg;
            
            %% Initialize the Variables
            obj.reluctance_torque=zeros(obj.size_gamma,obj.size_theta);
            
            %% Start a series of finite element caculations to Flux_PM of VPM
            openfemm(1);
            opendocument(obj.MyModel);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            mi_saveas('temp.fem');
            
            for gamma_deg = gamma_min:dgamma:gamma_max
                
                m = round((gamma_deg-gamma_min)/dgamma)+1;
                gamma_rad = gamma_deg*deg;
                MyIdCurrent = -obj.Im*sin(gamma_rad);     % direct current in phase current amplitude scaling
                MyIqCurrent = obj.Im*cos(gamma_rad);      % quadrature current phase current amplitude scaling
                
                for alpha_deg = alpha_min:dalpha:alpha_max
                    n = round((alpha_deg-alpha_min)/dalpha)+1;
                    % useful skills to caculate the time of circulation
                    starttime = clock;
                    
                    %% Set the initial rotor position
                    tta0 = initial_angle/obj.testmotor.rotor.pole_of_pairs_Rotor; %make the d axis is coinside with A axis
                    Rotor_angle=-(n-1)*dalpha+tta0; % rotor angle in degrees VPM clockwise rotation
                    mi_modifyboundprop('AGE',10,Rotor_angle);
                    % Set the Hc of PM material
                    mi_modifymaterial('PM',3,0);
                    
                    %% Control the current of phase
                    % make sure that the current is set to the appropriate value for this iteration.
                    tta = ((n-1)*dalpha)*deg*obj.testmotor.rotor.pole_of_pairs_Rotor;
                    Id = [cos(tta), cos(tta-2*pi/3), cos(tta+2*pi/3)];
                    Iq =-[sin(tta), sin(tta-2*pi/3), sin(tta+2*pi/3)];
                    Itot =  MyIdCurrent*Id + MyIqCurrent*Iq;
                    mi_setcurrent('A', Itot(1));
                    mi_setcurrent('B', Itot(2));
                    mi_setcurrent('C', Itot(3));
                    
                    %% Finite Elements Analysis Start
                    mi_analyze(1);
                    mi_loadsolution;
                    mo_smooth('off');  % flux smoothing algorithm is off to increase the speed of code.
                    
                    %% Magnetic Postprocessor Command
                    obj.reluctance_torque(m,n)=mo_gapintegral('AGE',0);  % Caculate the flux
                    
                    %% Close the Postprocessor and Delete the data of FEA
                    mo_close;
                    delete('temp.ans');
                    fprintf('% i of % i :: %f seconds ::  %f N*m \n',Rotor_angle*obj.MyLowestHarmonic*obj.testmotor.rotor.pole_of_pairs_Rotor,obj.theta_max_deg,etime(clock,starttime),obj.reluctance_torque(m,n));
                end
                fprintf('gamma = % i :: %f N*m \n', gamma_deg, mean(obj.reluctance_torque(m,:)));
            end
            
            %% Clean up after all finite element cacualation are finished
            closefemm;
            delete('temp.fem');
        end
        
        function Save_Data(obj)
            % Draw the 3D figure of the results
            [X,Y] = meshgrid(obj.theta_min_deg:obj.dtheta_deg:obj.theta_max_deg, obj.gamma_min_deg:obj.dgamma_deg:obj.gamma_max_deg);
            surf(X,Y,obj.reluctance_torque);
            xlabel('Advanced angle of current,Degrees');
            ylabel('Torque, N.m');
            title('Torque VS theta');
            grid on;
            % save the figure document and useful parameters
            figurename = strcat(obj.testmotor.name,'_Reluctance_Torque.fig');
            obj.figureaddress = strcat(obj.documentsaddress,'\', figurename);
            savefig(obj.figureaddress);        
            obj.testmotor.torque.reluctance_torque =obj.reluctance_torque;
        end
    end
end