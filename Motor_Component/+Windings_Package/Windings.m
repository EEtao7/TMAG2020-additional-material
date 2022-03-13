classdef Windings < handle
    % define the Windings class
    properties
        name
        precision_of_windings
        group_number_of_winding
        group_number_of_windingA
        group_number_of_windingB
        group_number_of_windingC
        pole_pairs_of_stator
        number_of_slots
        number_of_phase
        number_of_layers
          % define the stator where winding will be fixed
        stator
        % define the current of coils
        current_of_windings
    end
    properties(Constant)
        automesh = 1
        base_point = [0,0]
    end
    
    methods
        %          %% Constructor function
        %         function obj = Windings(name, precision_of_windings, group_number_of_windings, number_of_phase)
        %             obj.name = name;
        %             obj.precision_of_windings = precision_of_windings;
        %             obj.group_number_of_windings = group_number_of_windings;
        %             obj.number_of_phase = number_of_phase;
        %         end
    end
    %% Abstract methods
    methods(Abstract)
        define_material(obj);
        draw_single_element(obj);
    end
end