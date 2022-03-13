classdef AGE_Inner_Rotor < Air_gap_Package.AirGap
    % define the AGE_Inner_Rotor class
    properties   
        % define the airgap between stator and rotor
        stator
        AGE_flag = 6
        outer_AGE_angle = 0
        inner_AGE_angle = 0
    end
    
    %% define the Dependent properties
    properties(Dependent)
        outer_radius_of_AGE;
        inner_radius_of_AGE;
    end
    
    methods
%         %% Constructor function
%         function obj = AGE_Inner_Rotor(name,group_number_of_AirGap,precision_of_AirGap,air_gap_length, stator, AGE_flag, outer_AGE_angle, inner_AGE_angle)
%             obj = obj@AirGap(name,group_number_of_AirGap,precision_of_AirGap,air_gap_length);
%             obj.stator = stator;
%             obj.AGE_flag = AGE_flag;
%             obj.outer_AGE_angle = outer_AGE_angle;
%             obj.inner_AGE_angle = inner_AGE_angle;
%         end

        %% Dependent properties caculation function
        function outer_radius_of_AGE = get.outer_radius_of_AGE(obj)
            outer_radius_of_AGE = obj.stator.inner_radius-0.25*obj.air_gap_length;
        end
        
        function inner_radius_of_AGE = get.inner_radius_of_AGE(obj)
            inner_radius_of_AGE =  obj.stator.inner_radius-0.75*obj.air_gap_length;
        end
        
        %% Various function of femm
        function define_material(obj)
            % Properties of Material
                % define the properties of Air
                mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
            % Properties of Boundary
                % define the AGE Boundary Condition
                mi_addboundprop('AGE',0, 0, 0, 0, 0, 0, 0, 0, obj.AGE_flag, obj.inner_AGE_angle, obj.outer_AGE_angle);
        end
        
        function draw_single_element(obj)
            %% Draw the AGE
            mi_drawarc(point(obj.outer_radius_of_AGE,0),point(obj.outer_radius_of_AGE,180),180,obj.precision_of_AirGap);
            mi_selectarcsegment(point(obj.outer_radius_of_AGE, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_AirGap,'AGE',0,obj.group_number_of_AirGap);
            mi_clearselected;
            mi_drawarc(point(obj.outer_radius_of_AGE,180),point(obj.outer_radius_of_AGE,0),180,obj.precision_of_AirGap);
            mi_selectarcsegment(point(obj.outer_radius_of_AGE, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_AirGap,'AGE',0,obj.group_number_of_AirGap);
            mi_clearselected;
            
            mi_drawarc(point(obj.inner_radius_of_AGE,0),point(obj.inner_radius_of_AGE,180),180,obj.precision_of_AirGap);
            mi_selectarcsegment(point(obj.inner_radius_of_AGE, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_AirGap,'AGE',0,obj.group_number_of_AirGap);
            mi_clearselected;
            mi_drawarc(point(obj.inner_radius_of_AGE,180),point(obj.inner_radius_of_AGE,0),180,obj.precision_of_AirGap);
            mi_selectarcsegment(point(obj.inner_radius_of_AGE, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_AirGap,'AGE',0,obj.group_number_of_AirGap);
            mi_clearselected;
            
            %% Set the material of air-gap & AGE
            inner_air_point = point(obj.inner_radius_of_AGE-0.1*obj.air_gap_length,0);
            mi_addblocklabel(inner_air_point);
            mi_selectlabel(inner_air_point);
            mi_setblockprop('Air', obj.automesh, obj.precision_of_AirGap, '<None>', 0, obj.group_number_of_AirGap, 0);
            mi_clearselected;
            outer_air_point = point(obj.outer_radius_of_AGE+0.1*obj.air_gap_length,0);
            mi_addblocklabel(outer_air_point);
            mi_selectlabel(outer_air_point);
            mi_setblockprop('Air', obj.automesh, obj.precision_of_AirGap, '<None>', 0, obj.group_number_of_AirGap, 0);
            mi_clearselected;
            AGE_point = point(0.5*(obj.inner_radius_of_AGE+obj.outer_radius_of_AGE),0);
            mi_addblocklabel(AGE_point);
            mi_selectlabel(AGE_point);
            mi_setblockprop('<No Mesh>', obj.automesh, obj.precision_of_AirGap, '<None>', 0, obj.group_number_of_AirGap,0);
            mi_clearselected;
           
        end    
    end
end