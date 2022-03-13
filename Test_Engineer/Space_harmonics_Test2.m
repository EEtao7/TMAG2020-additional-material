classdef Space_harmonics_Test2 < TestProject
    properties
        % test items(测试项目)
        MyModel
        airgap = 'AGE'
        B_air_gap
        B_air_gap_fft
        documentsaddress = 'E:\My_Study\ZJU_master_study\master research\Low speed servo motor\Matlab_and_Simulink_documents\figure documents\Space harmonics';
        figureaddress
        
        % adjust test parameters（可调节的测试参数）
        initial_angle
        alpha_min_deg
        alpha_max_deg
        dalpha_deg
        
    end
    properties(Dependent)
        size_alpha
    end
    methods
        function size_alpha = get.size_alpha(obj)
            size_alpha = round((obj.alpha_max_deg-obj.alpha_min_deg)/obj.dalpha_deg)+1;
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
            %             Periodic_coeff = 1; % when caculating the Inductance and Flux by femm Periodic BC, you need to multiplicate it. Not to the torque in femm.
            %             Park_coeff = 1.5; % when caculating the torque and power by equations in dq coordinate system.
            
            %% Other User-defined  Parameters
            % Set the mechanical angle alpha
            alpha_min = obj.alpha_min_deg;
            alpha_max = obj.alpha_max_deg;
            dalpha = obj.dalpha_deg;
            
            %% Initialize the Variables
            obj.B_air_gap=zeros(1,obj.size_alpha);
            
            %% Start a series of finite element caculations to Flux_PM of VPM
            % useful skills to caculate the time of circulation
            starttime = clock;
            openfemm(1);
            opendocument(obj.MyModel);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            mi_saveas('temp.fem');
            
            % Set the current in the circuit coils
            mi_setcurrent('A',0);
            mi_setcurrent('B',0);
            mi_setcurrent('C',0);
            
            %% Set the initial rotor position
            tta0 = obj.initial_angle; %make the d axis is coinside with A axis
            
            Rotor_angle = tta0; % rotor angle in degrees VPM clockwise rotation
            mi_modifyboundprop('AGE',10,Rotor_angle);% inner rotor
            Stator_angle = 360/obj.testmotor.stator.number_of_slots/2;
            mi_modifyboundprop('AGE',11,Stator_angle);% inner rotor
            
            %% Finite Elements Analysis Start
            mi_analyze(1);
            mi_loadsolution;
            mo_smooth('off');  % flux smoothing algorithm is off to increase the speed of code.
            
            %% Magnetic Postprocessor Command
            %% Recoed the B in VPM air-gap
            
            for alpha = alpha_min : dalpha : alpha_max
                m = round((alpha-alpha_min)/dalpha+1);
                B = mo_getgapb(obj.airgap, alpha);
                obj.B_air_gap(1,m)=B(1,1);
            end
            
            %% Close the Postprocessor and Delete the data of FEA
            mo_close;
            delete('temp.ans');
            fprintf('complete the magnetic field caculation :: %f seconds \n', etime(clock,starttime));
            
            %% Clean up after all finite element cacualation are finished
            closefemm;
            delete('temp.fem');
            
            
        end
        
        function Save_Data(obj)
            % Draw the figure of the results
            subplot(2,1,1)
            plot(obj.alpha_min_deg:obj.dalpha_deg:obj.alpha_max_deg, obj.B_air_gap);
            xlabel('Mechanical Angle, Degrees');
            ylabel('Amplitude of Br,Tesla');
            title('Br of air-gap in VPM');
            grid on;
            
            % Compute the square of the amplitude of each harmonic at the centroid of
            % each element in the mesh. Matlab's built-in FFT function makes this easy.
            subplot(2,1,2)
            sample_frequency = obj.size_alpha;
            [f, obj.B_air_gap_fft]=FFT(obj.B_air_gap, sample_frequency);
            bar(f, obj.B_air_gap_fft);
            xlabel('Harmonics Order');
            ylabel('Amplitude');
            title('Br of air-gap in VPM (FFT)');
            grid on;
            
            % save the figure document and useful parameters
            figurename = strcat(obj.testmotor.name,'_Space_harmonics.fig');
            obj.figureaddress = strcat(obj.documentsaddress,'\', figurename);
            savefig(obj.figureaddress);
        end
    end
end