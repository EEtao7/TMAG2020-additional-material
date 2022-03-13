classdef MotorDesignEngineer < handle
    % MotorDesign who can use femm to construct builder.
    properties
        builder     % define the builder which will be constructed by MotorDesignEngineer.
    end
    methods
        function construct(obj)
            % define the sequence of construct a builder
            obj.builder.define_material();
            obj.builder.draw_single_element();                
        end
    end  
end