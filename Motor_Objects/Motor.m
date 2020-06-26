classdef Motor < handle
    % define the Motor class to be design, test and control
    properties 
        documentsaddress = 'E:\My_Study\ZJU_master_study\master research\Low speed servo motor\Matlab_and_Simulink_documents\femm documents';
        name
        fileaddress
        stator
        air_gap
        rotor
        windings
        stack_length
        controlparameters = ControlParameters();
        torque = Torque();
    end
    
    methods
        function create_femm(obj)
            % open femm and save as a femm document
            openfemm(1);
            newdocument(0);
            mi_smartmesh(0);  % use relatively coarse mesh to save time
            
            % define the name of femm document
            obj.name = strcat(obj.stator.name, obj.rotor.name,...
                num2str(obj.stator.pole_pairs_of_stator),'-', num2str(obj.rotor.pole_of_pairs_Rotor),'-', num2str(obj.stator.number_of_slots*obj.stator.number_of_auxiliary_tooth));
%             obj.name = 'TestMotor';  % ²âÊÔÄ£Ê½

            file = strcat(obj.name,'.fem');
            % define the file address of femm document 
            obj.fileaddress = strcat(obj.documentsaddress,'\', file);
            mi_saveas(obj.fileaddress);
            
            % define the parameters of problem solve
            mi_probdef(0, 'millimeters', 'planar', 1e-8, obj.stack_length, 30, 0);
        end
        
        function close_femm(obj)
            % save the femm document and close the femm
            mi_saveas(obj.fileaddress);
            closefemm; 
            % save the parameters of the motor
            
        end
        
    end
end