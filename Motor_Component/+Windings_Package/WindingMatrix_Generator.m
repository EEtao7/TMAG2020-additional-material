classdef WindingMatrix_Generator < handle
    % define the Windings class
    properties
        name
        number_of_phase
        number_of_slots    
        number_of_layers
        pole_pairs_of_windings
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
        define_name(obj);
        generate_the_windingmatrix_phaseA(obj);
    end
end