classdef AirGap < handle
    % define the AirGap class
    properties
        name
        group_number_of_AirGap
        precision_of_AirGap 
        air_gap_length
    end
    properties(Constant)
        automesh = 1;
        base_point = [0,0];
    end

    methods
%          %% Constructor function
%         function obj = AirGap(name,group_number_of_AirGap,precision_of_AirGap,air_gap_length)
%             obj.name = name;
%             obj.group_number_of_AirGap = group_number_of_AirGap;
%             obj.precision_of_AirGap = precision_of_AirGap;
%             obj.air_gap_length = air_gap_length;
%         end
    end
    %% Abstract methods
    methods(Abstract)
         define_material(obj);
        draw_single_element(obj);
    end
end