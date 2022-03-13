classdef OuterStator < handle
    % define the Stator class
    properties
        name
        group_number_of_stator
        precision_of_stator_iron
        outer_radius
    end
    properties(Constant)
        automesh = 1
        base_point = [0,0]
    end
    
    methods
        %          %% Constructor function
        %         function obj = Stator(name, group_number_of_stator, precision_of_stator_iron, outer_radious)
        %             obj.name = name;
        %             obj.group_number_of_stator = group_number_of_stator;
        %             obj.precision_of_stator_iron = precision_of_stator_iron;
        %             obj.outer_radious = outer_radious;
        %         end
    end
    %% Abstract methods
    methods(Abstract)
        define_material(obj);
        draw_single_element(obj);
    end
end