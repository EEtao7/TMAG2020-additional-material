classdef MagneticFigureGenerator < TestProject
    properties
        MyModel
        % adjust parameters
        Motor_flag
        Rotor_flag
        Figure_flag
        Im
        initial_angle_deg
        
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
            Pi=pi;
            deg=Pi/180.;
            
            % some useful coefficients to transfer the values
            % Periodic_coeff = 1; % when caculating the Inductance and Flux by femm Periodic BC, you need to multiplicate it. Not to the torque in femm.
            %             Park_coeff = 1.5; % when caculating the torque and power by equations in dq coordinate system.
            
            %% Other User-defined  Parameters
            % Set the mechanical angle alpha
            initial_angle = 0;
            alpha_min = obj.theta_min_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            alpha_max = obj.theta_max_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            dalpha = obj.dtheta_deg/obj.testmotor.rotor.pole_of_pairs_Rotor;
            
            
            %% Start a series of finite element caculations for the MagneticFigure ;
            openfemm(0);
            opendocument(obj.MyModel);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            mi_saveas('temp.fem');
            
            if(obj.Figure_flag == 1)
                %% MagneticFigure of Noload
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
                    
                    % save the magneticfield figure
                    mo_showdensityplot(0, 0, 3, 0, 'mag');
                    x1 = -105;
                    y1 = -105;
                    x2 = 105;
                    y2 = 105;
                    address = 'E:\My_Picture\';
                    name = 'test Motor_Noload';
                    number = m;
                    MagneticfieldFigure(x1, y1, x2, y2, address, name, number);
                    %% Close the Postprocessor and Delete the data of FEA
                    mo_close;
                    delete('temp.ans');
                    fprintf('% i of % i :: %f seconds, % i \n',Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor,obj.theta_max_deg,etime(clock,starttime),obj.Figure_flag);
                end
                
            elseif(obj.Figure_flag == 2)
                %% MagneticFigure of Load
                % Set the current in the circuit coils
                MyIdCurrent = -obj.Im*sin(0);     % direct current in phase current amplitude scaling
                MyIqCurrent = obj.Im*cos(0);      % quadrature current phase current amplitude scaling
                
                for alpha_deg = alpha_min:dalpha:alpha_max
                    m = round((alpha_deg-alpha_min)/dalpha)+1;
                    % useful skills to caculate the time of circulation
                    starttime = clock;
                    
                    %% Set the initial rotor position
                    tta0 = obj.initial_angle_deg/obj.testmotor.rotor.pole_of_pairs_Rotor; %make the d axis is coinside with A axis
                    if (obj.Motor_flag==1)
                        Rotor_angle=-(m-1)*dalpha-tta0; % rotor angle in degrees  clockwise rotation
                    elseif(obj.Motor_flag==2)
                        Rotor_angle=(m-1)*dalpha+tta0; % rotor angle in degrees  anticlockwise rotation
                    end
                    if(obj.Rotor_flag==1)
                        mi_modifyboundprop('AGE',11,Rotor_angle);% outer rotor
                    elseif(obj.Rotor_flag==2)
                        mi_modifyboundprop('AGE',10,Rotor_angle);% inner rotor
                    end
                    
                    %% Control the current of phase
                    % make sure that the current is set to the appropriate value for this iteration.
                    tta = ((m-1)*dalpha)*deg*obj.testmotor.rotor.pole_of_pairs_Rotor;
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
                    
                    % save the magneticfield figure
                    mo_showdensityplot(0, 0, 3, 0, 'mag');
                    x1 = -105;
                    y1 = -105;
                    x2 = 105;
                    y2 = 105;
                    address = 'E:\My_Picture\';
                    name = 'test Motor_Load';
                    number = m;
                    MagneticfieldFigure(x1, y1, x2, y2, address, name, number);
                    %% Close the Postprocessor and Delete the data of FEA
                    mo_close;
                    delete('temp.ans');
                    fprintf('% i of % i :: %f seconds, % i \n',Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor,obj.theta_max_deg,etime(clock,starttime),obj.Figure_flag);
                end
                
            elseif(obj.Figure_flag == 3)
                %% MagneticFigure of NoPM
                % Set the Hc of PM material
                mi_modifymaterial('PM',3,0);
                MyIdCurrent = -obj.Im*sin(0);     % direct current in phase current amplitude scaling
                MyIqCurrent = obj.Im*cos(0);      % quadrature current phase current amplitude scaling
                
                for alpha_deg = alpha_min:dalpha:alpha_max
                    m = round((alpha_deg-alpha_min)/dalpha)+1;
                    % useful skills to caculate the time of circulation
                    starttime = clock;
                    
                    %% Set the initial rotor position
                    tta0 = obj.initial_angle_deg/obj.testmotor.rotor.pole_of_pairs_Rotor; %make the d axis is coinside with A axis
                    if (obj.Motor_flag==1)
                        Rotor_angle=-(m-1)*dalpha-tta0; % rotor angle in degrees  clockwise rotation
                    elseif(obj.Motor_flag==2)
                        Rotor_angle=(m-1)*dalpha+tta0; % rotor angle in degrees  anticlockwise rotation
                    end
                    if(obj.Rotor_flag==1)
                        mi_modifyboundprop('AGE',11,Rotor_angle);% outer rotor
                    elseif(obj.Rotor_flag==2)
                        mi_modifyboundprop('AGE',10,Rotor_angle);% inner rotor
                    end
                    
                    %% Control the current of phase
                    % make sure that the current is set to the appropriate value for this iteration.
                    tta = ((m-1)*dalpha)*deg*obj.testmotor.rotor.pole_of_pairs_Rotor;
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
                    
                    % save the magneticfield figure
                    mo_showdensityplot(0, 0, 3, 0, 'mag');
                    x1 = -105;
                    y1 = -105;
                    x2 = 105;
                    y2 = 105;
                    address = 'E:\My_Picture\';
                    name = 'test Motor_NoPM';
                    number = m;
                    MagneticfieldFigure(x1, y1, x2, y2, address, name, number);
                    %% Close the Postprocessor and Delete the data of FEA
                    mo_close;
                    delete('temp.ans');
                    fprintf('% i of % i :: %f seconds, % i \n',Rotor_angle*obj.testmotor.rotor.pole_of_pairs_Rotor,obj.theta_max_deg,etime(clock,starttime),obj.Figure_flag);
                end
            end
            %% Clean up after all finite element cacualation are finished
            closefemm;
            delete('temp.fem');
        end
        function Save_Data(obj)
            fprintf('Figure % i Mission Complete\n',obj.Figure_flag);
        end
    end
end