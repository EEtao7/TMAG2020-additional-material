classdef Inner_Rotor < handle
    % define the Inner Rotor class
    properties
        name
        group_number_of_rotor
        precision_of_rotor_iron
        inner_radious
    end
    properties(Constant)
        automesh = 1
        base_point = [0,0]
    end
    
    methods
        %          %% Constructor function
        %         function obj = Rotor(name,number_of_group,precision,inner_radious)
        %             obj.name = name;
        %             obj.group_number_of_rotor = number_of_group;
        %             obj.precision_of_rotor_iron = precision;
        %             obj.inner_radious = inner_radious;
        %         end
    end
    %% Abstract methods
    methods(Abstract)
        define_material(obj);
        draw_single_element(obj);
    end
end