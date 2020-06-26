classdef SurfacemountedPM_Inner_Rotor < Rotor_Package.Inner_Rotor
    % define the Surface mounted PM Rotor class
    properties
        pole_of_pairs_Rotor
        thickness_of_PM
        radius_of_PM
        ratio_of_rotor_PM_bottom
        ratio_of_rotor_PM_top
        precision_of_PM
        group_number_of_PM
        material_of_PM
        material_of_Iron
        nonlinear_flag
    end
    
    %% define the Dependent properties
    properties(Dependent)
        angle_region_of_PM;
        angle_arc_of_PM_bottom;
        angle_arc_of_PM_top
        outer_radious;
    end
    
    methods
        %         %% Constructor function
        %         function obj = SurfacemountedPM(name, group_number_of_rotor, precision_of_rotor_iron,  inner_radious, pole_of_pairs_Rotor, thick_of_rotor_PM, radius_of_PM, ratio_of_rotor_PM,precision_of_PM,group_number_of_PM, material_of_PM, material_of_Iron)
        %             obj = obj@Rotor(name, group_number_of_rotor, precision_of_rotor_iron, inner_radious);
        %             obj.pole_of_pairs_Rotor = pole_of_pairs_Rotor;
        %             obj.thick_of_rotor_PM = thick_of_rotor_PM;
        %             obj.radius_of_PM = radius_of_PM;
        %             obj.ratio_of_rotor_PM = ratio_of_rotor_PM;
        %             obj.precision_of_PM = precision_of_PM;
        %             obj.group_number_of_PM = group_number_of_PM;
        %             obj.material_of_PM = material_of_PM;
        %             obj.material_of_Iron = material_of_Iron;
        %         end
        
        %% Dependent properties caculation function
        function angle_region_of_PM = get.angle_region_of_PM(obj)
            angle_region_of_PM = 360/(obj.pole_of_pairs_Rotor*2);
        end
        
        function angle_arc_of_PM_bottom = get.angle_arc_of_PM_bottom(obj)
            angle_arc_of_PM_bottom = obj.ratio_of_rotor_PM_bottom*obj.angle_region_of_PM;
        end
        
         function angle_arc_of_PM_top = get.angle_arc_of_PM_top(obj)
            angle_arc_of_PM_top = obj.ratio_of_rotor_PM_top*obj.angle_region_of_PM;
        end
        
        function outer_radious = get.outer_radious(obj)
            outer_radious = obj.radius_of_PM-obj.thickness_of_PM;
        end
        
        %% Various function of femm
        function define_material(obj)
            % Properties of Material
            % define the properties of Air
            mi_addmaterial('Axis', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
            % define the properties of PM
            mi_addmaterial('PM',1.05, 1.05, obj.material_of_PM.PM_Hc, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
            % define the properties of Iron
            mi_addmaterial('Rotor_Iron',1000000, 1000000, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
            if(obj.nonlinear_flag)
                mi_addbhpoints('Rotor_Iron', obj.material_of_Iron.bhcurve);
            end
            % Properties of Boundary
            % define the A = 0 Boundary Condition
            mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);
        end
        
        function draw_single_element(obj)
            %% surface-mounted PM of rotor
            % draw the axis and rotor of VPM
            mi_drawarc(point(obj.inner_radious,0),point(obj.inner_radious,180),180,obj.precision_of_rotor_iron);
            mi_selectarcsegment(point(obj.inner_radious, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_rotor_iron,'A=0',1,obj.group_number_of_rotor);
            mi_clearselected;
            mi_drawarc(point(obj.inner_radious,180),point(obj.inner_radious,0),180,obj.precision_of_rotor_iron);
            mi_selectarcsegment(point(obj.inner_radious, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_rotor_iron,'A=0',1,obj.group_number_of_rotor);
            mi_clearselected;
            
            mi_drawarc(point(obj.outer_radious,0),point(obj.outer_radious,180),180,obj.precision_of_rotor_iron);
            mi_selectarcsegment(point(obj.outer_radious, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_rotor_iron,'<None>',0,obj.group_number_of_rotor);
            mi_clearselected;
            mi_drawarc(point(obj.outer_radious,180),point(obj.outer_radious,0),180,obj.precision_of_rotor_iron);
            mi_selectarcsegment(point(obj.outer_radious, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_rotor_iron,'<None>',0,obj.group_number_of_rotor);
            mi_clearselected;
            
            % draw the geometry of rotor PM
            mi_drawline(point(obj.outer_radious,-obj.angle_arc_of_PM_bottom*0.5),point(obj.radius_of_PM,-obj.angle_arc_of_PM_top*0.5));
            mi_drawline(point(obj.outer_radious,obj.angle_arc_of_PM_bottom*0.5),point(obj.radius_of_PM,obj.angle_arc_of_PM_top*0.5));
            mi_drawarc(point(obj.radius_of_PM,-obj.angle_arc_of_PM_top*0.5),point(obj.radius_of_PM,obj.angle_arc_of_PM_top*0.5),obj.angle_arc_of_PM_top,obj.precision_of_PM);
            
            % select and set properties of geometry of rotor PM
            mi_selectsegment(0.5*(point(obj.outer_radious,-obj.angle_arc_of_PM_bottom*0.5)+point(obj.radius_of_PM,-obj.angle_arc_of_PM_top*0.5)));
            mi_selectsegment(0.5*(point(obj.outer_radious,obj.angle_arc_of_PM_bottom*0.5)+point(obj.radius_of_PM,obj.angle_arc_of_PM_top*0.5)));
            mi_setsegmentprop('<None>',obj.precision_of_PM,obj.automesh,0,obj.group_number_of_PM);
            mi_clearselected;
            mi_selectarcsegment(point(obj.radius_of_PM,0.5*((obj.angle_arc_of_PM_top*0.5)+(-obj.angle_arc_of_PM_top*0.5))));
            mi_setarcsegmentprop(obj.precision_of_PM,'<None>',0,obj.group_number_of_PM);
            mi_clearselected;
            % copy the PM of rotor
            mi_selectgroup(1);
            mi_copyrotate(obj.base_point,obj.angle_region_of_PM,(2*obj.pole_of_pairs_Rotor-1));
            mi_clearselected;
            
            %% Set and Define the BlockLabel
            % set the material of axis
            axis_point = obj.base_point;
            mi_addblocklabel(axis_point);
            mi_selectlabel(axis_point);
            mi_setblockprop('Axis',obj.automesh,obj.precision_of_rotor_iron,'<None>',0,obj.group_number_of_rotor,0);
            mi_clearselected;
            
            % set the material of Iron
            rotor_Iron_point = point(0.5*(obj.inner_radious+obj.outer_radious),0);
            mi_addblocklabel(rotor_Iron_point);
            mi_selectlabel(rotor_Iron_point);
            mi_setblockprop('Rotor_Iron',obj.automesh,obj.precision_of_rotor_iron,'<None>',0,obj.group_number_of_rotor,0);
            mi_clearselected;
            
            % surface-mounted PM of rotor
            PM_N_point = point((obj.outer_radious+0.5*obj.thickness_of_PM),0);
            mi_addblocklabel(PM_N_point);
            mi_selectlabel(PM_N_point);
            mi_setblockprop('PM',obj.automesh,obj.precision_of_PM,'<None>',0,obj.group_number_of_PM,0);
            mi_copyrotate(obj.base_point,2*obj.angle_region_of_PM,obj.pole_of_pairs_Rotor-1);
            mi_clearselected;
            PM_S_point = point((obj.outer_radious+0.5*obj.thickness_of_PM),obj.angle_region_of_PM);
            mi_addblocklabel(PM_S_point);
            mi_selectlabel(PM_S_point);
            mi_setblockprop('PM',obj.automesh,obj.precision_of_PM,'<None>',180+obj.angle_region_of_PM,obj.group_number_of_PM,0);
            mi_copyrotate(obj.base_point,2*obj.angle_region_of_PM,obj.pole_of_pairs_Rotor-1);
            mi_clearselected;
            
        end
    end
end