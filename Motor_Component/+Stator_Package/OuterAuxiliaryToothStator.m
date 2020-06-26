classdef OuterAuxiliaryToothStator < Stator_Package.OuterStator
    % define the outer auxiliary tooth stator class
    properties
        % poles parameters of stator
        pole_pairs_of_stator
        number_of_slots
        number_of_auxiliary_tooth
        % tangential parameters of stator
        ratio_of_tooth
        ratio_of_auxiliary_tooth
        % radial parameters of stator
        slot_depth
        auxiliary_slot_depth
        thick_of_yoke_stator
        thick_of_auxiliary_yoke_stator
        adjust_coefficient_of_thickness_of_auxiliary_yoke_stator
        % parameters of stator
        material_of_Iron
        nonlinear_flag
    end
    
    %% define the Dependent properties
    properties(Dependent)
        number_of_modulation_rings
        % angle region of stator
        angle_region_of_tooth
        angle_arc_of_tooth
        angle_region_of_auxiliary_tooth
        angle_arc_of_auxiliary_tooth
        angle_arc_of_auxiliary_slot
        % radial parameters of stator
        outer_radius_of_slot
        inner_radius_of_slot
        adjust_radius_of_slot
        radius_of_auxiliary_slot
        inner_radius
    end
    
    methods
        %         %% Constructor function
        %         function obj = AuxiliaryToothStator(name, group_number_of_stator, precision_of_stator_iron, outer_radius, pole_of_pairs_stator, number_of_slots, number_of_auxiliary_tooth, ratio_of_tooth, ratio_of_auxiliary_tooth, slot_depth, auxiliary_slot_depth, thick_of_yoke_stator, thick_of_auxiliary_yoke_stator, material_of_Iron)
        %             obj = obj@Stator(name, group_number_of_stator, precision_of_stator_iron, outer_radius);
        %             obj.pole_of_pairs_stator = pole_of_pairs_stator;
        %             obj.number_of_slots = number_of_slots;
        %             obj.number_of_auxiliary_tooth = number_of_auxiliary_tooth;
        %             obj.ratio_of_tooth = ratio_of_tooth;
        %             obj.ratio_of_auxiliary_tooth = ratio_of_auxiliary_tooth;
        %             obj.slot_depth = slot_depth;
        %             obj.auxiliary_slot_depth = auxiliary_slot_depth;
        %             obj.thick_of_yoke_stator = thick_of_yoke_stator;
        %             obj.thick_of_auxiliary_yoke_stator = thick_of_auxiliary_yoke_stator;
        %             obj.material_of_Iron = material_of_Iron;
        %         end
        
        %% Dependent properties caculation function
        function number_of_modulation_rings = get.number_of_modulation_rings(obj)
            number_of_modulation_rings = obj.number_of_slots*obj.number_of_auxiliary_tooth;
        end
        function angle_region_of_tooth = get.angle_region_of_tooth(obj)
            angle_region_of_tooth = 360/obj.number_of_slots;
        end
        function angle_arc_of_tooth = get.angle_arc_of_tooth(obj)
            angle_arc_of_tooth = obj.ratio_of_tooth*obj.angle_region_of_tooth;
        end
        function angle_region_of_auxiliary_tooth = get.angle_region_of_auxiliary_tooth(obj)
            angle_region_of_auxiliary_tooth = obj.angle_region_of_tooth/obj.number_of_auxiliary_tooth;
        end
        function angle_arc_of_auxiliary_tooth = get.angle_arc_of_auxiliary_tooth(obj)
            angle_arc_of_auxiliary_tooth = obj.ratio_of_auxiliary_tooth*obj.angle_region_of_auxiliary_tooth;
        end
        function angle_arc_of_auxiliary_slot = get.angle_arc_of_auxiliary_slot(obj)
            angle_arc_of_auxiliary_slot = (1-obj.ratio_of_auxiliary_tooth)*obj.angle_region_of_auxiliary_tooth;
        end
        function outer_radius_of_slot = get.outer_radius_of_slot(obj)
            outer_radius_of_slot = obj.outer_radius-obj.thick_of_yoke_stator;
        end
        function inner_radius_of_slot = get.inner_radius_of_slot(obj)
            inner_radius_of_slot = obj.outer_radius_of_slot-obj.slot_depth;
        end
        function adjust_radius_of_slot = get.adjust_radius_of_slot(obj)
            adjust_radius_of_slot = obj.outer_radius_of_slot-obj.slot_depth-obj.thick_of_auxiliary_yoke_stator*(1-obj.adjust_coefficient_of_thickness_of_auxiliary_yoke_stator);
        end
        function radius_of_auxiliary_slot = get.radius_of_auxiliary_slot(obj)
            radius_of_auxiliary_slot = obj.inner_radius_of_slot-obj.thick_of_auxiliary_yoke_stator;
        end
        function inner_radius = get.inner_radius(obj)
            inner_radius = obj.radius_of_auxiliary_slot-obj.auxiliary_slot_depth;
        end
        
        %% Various function of femm
        function define_material(obj)
            % Properties of Material
            % define the properties of Iron
            mi_addmaterial('Stator_Iron',1000000, 1000000, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
            if (obj.nonlinear_flag)
                mi_addbhpoints('Stator_Iron', obj.material_of_Iron.bhcurve);
            end
            % Properties of Boundary
            % define the A = 0 Boundary Condition
            mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);
            % define the Periodic Boundary Condition
            mi_addboundprop('Periodic BC1', 0, 0, 0, 0, 0, 0, 0, 0, 4);
            mi_addboundprop('Periodic BC2', 0, 0, 0, 0, 0, 0, 0, 0, 4);
            % define the Point Properties
            mi_addpointprop('A=0', 0,0);
        end
        
        function draw_single_element(obj)
            %% Draw the Staor of VPM
            % auxiliary slots of split tooth stator
            %% draw the geometry of stator
            % draw the auxiliary tooth
            if(mod(obj.number_of_auxiliary_tooth,2) == 1)
                % Na is odd
                base_angle = obj.angle_arc_of_auxiliary_tooth*0.5;
                mi_drawarc(point(obj.inner_radius,0),point(obj.inner_radius,base_angle),obj.angle_arc_of_auxiliary_tooth*0.5,0.2*obj.precision_of_stator_iron);
                mi_selectarcsegment(point(obj.inner_radius,0.5*(0+(base_angle))));
                mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
                mi_clearselected;
                
                mi_drawline(point(obj.inner_radius,base_angle),point(obj.radius_of_auxiliary_slot,base_angle));
                mi_selectsegment(0.5*(point(obj.inner_radius,base_angle)+point(obj.radius_of_auxiliary_slot,base_angle)));
                mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                mi_clearselected;
                
                if (obj.number_of_auxiliary_tooth > 1)
                    for m = 1:fix(obj.number_of_auxiliary_tooth/2)
                        mi_drawline(point(obj.radius_of_auxiliary_slot,base_angle+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.radius_of_auxiliary_slot,base_angle+(m-1)*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawline(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawarc(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.inner_radius,base_angle+obj.angle_region_of_auxiliary_tooth+(m-1)*obj.angle_region_of_auxiliary_tooth),obj.angle_arc_of_auxiliary_tooth,0.2*obj.precision_of_stator_iron)
                        mi_selectarcsegment(point(obj.inner_radius,0.5*((base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)+(base_angle+obj.angle_region_of_auxiliary_tooth+(m-1)*obj.angle_region_of_auxiliary_tooth))));
                        mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawline(point(obj.inner_radius,base_angle+m*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+m*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.inner_radius,base_angle+m*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+m*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                    end
                end
                
            elseif(mod(obj.number_of_auxiliary_tooth,2) == 0)
                % Na is even
                base_angle = obj.angle_arc_of_auxiliary_slot*0.5+obj.angle_arc_of_auxiliary_tooth;
                mi_drawline(point(obj.radius_of_auxiliary_slot,0),point(obj.radius_of_auxiliary_slot,obj.angle_arc_of_auxiliary_slot*0.5));
                mi_selectsegment(0.5*(point(obj.radius_of_auxiliary_slot,0)+point(obj.radius_of_auxiliary_slot,obj.angle_arc_of_auxiliary_slot*0.5)));
                mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                mi_clearselected;
                
                mi_drawline(point(obj.inner_radius,obj.angle_arc_of_auxiliary_slot*0.5),point(obj.radius_of_auxiliary_slot,obj.angle_arc_of_auxiliary_slot*0.5));
                mi_selectsegment(0.5*(point(obj.inner_radius,obj.angle_arc_of_auxiliary_slot*0.5)+point(obj.radius_of_auxiliary_slot,obj.angle_arc_of_auxiliary_slot*0.5)));
                mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                mi_clearselected;
                
                mi_drawarc(point(obj.inner_radius,obj.angle_arc_of_auxiliary_slot*0.5),point(obj.inner_radius,base_angle),obj.angle_arc_of_auxiliary_tooth,0.2*obj.precision_of_stator_iron);
                mi_selectarcsegment(point(obj.inner_radius,0.5*((obj.angle_arc_of_auxiliary_slot*0.5)+(base_angle))));
                mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
                mi_clearselected;
                
                mi_drawline(point(obj.inner_radius,base_angle),point(obj.radius_of_auxiliary_slot,base_angle));
                mi_selectsegment(0.5*(point(obj.inner_radius,base_angle)+point(obj.radius_of_auxiliary_slot,base_angle)));
                mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                mi_clearselected;
                
                if (obj.number_of_auxiliary_tooth > 2)
                    for m = 1:fix(obj.number_of_auxiliary_tooth/2-1)
                        mi_drawline(point(obj.radius_of_auxiliary_slot,base_angle+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.radius_of_auxiliary_slot,base_angle+(m-1)*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawline(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawarc(point(obj.inner_radius,base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth),point(obj.inner_radius,base_angle+obj.angle_region_of_auxiliary_tooth+(m-1)*obj.angle_region_of_auxiliary_tooth),obj.angle_arc_of_auxiliary_tooth,0.2*obj.precision_of_stator_iron)
                        mi_selectarcsegment(point(obj.inner_radius,0.5*((base_angle+obj.angle_arc_of_auxiliary_slot+(m-1)*obj.angle_region_of_auxiliary_tooth)+(base_angle+obj.angle_region_of_auxiliary_tooth+(m-1)*obj.angle_region_of_auxiliary_tooth))));
                        mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
                        mi_clearselected;
                        
                        mi_drawline(point(obj.inner_radius,base_angle+m*obj.angle_region_of_auxiliary_tooth),point(obj.radius_of_auxiliary_slot,base_angle+m*obj.angle_region_of_auxiliary_tooth));
                        mi_selectsegment(0.5*(point(obj.inner_radius,base_angle+m*obj.angle_region_of_auxiliary_tooth)+point(obj.radius_of_auxiliary_slot,base_angle+m*obj.angle_region_of_auxiliary_tooth)));
                        mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
                        mi_clearselected;
                    end
                end
            end
            % draw the  tooth
            mi_drawline(point(obj.radius_of_auxiliary_slot,0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot)),point(obj.adjust_radius_of_slot,0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot)));
            mi_selectsegment(0.5*(point(obj.radius_of_auxiliary_slot,0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot))+point(obj.inner_radius_of_slot,0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot))));
            mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawline(point(obj.outer_radius_of_slot,0.5*obj.angle_arc_of_tooth)-[obj.slot_depth, 0], point(obj.adjust_radius_of_slot,0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot)));
            mi_selectsegment(point(obj.inner_radius_of_slot,0.5*((0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot))+(0.5*obj.angle_arc_of_tooth))));
            mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawline(point(obj.outer_radius_of_slot,0.5*obj.angle_arc_of_tooth)-[obj.slot_depth, 0], point(obj.outer_radius_of_slot,0.5*obj.angle_arc_of_tooth));
            mi_selectsegment(0.5*(point(obj.inner_radius_of_slot,0.5*obj.angle_arc_of_tooth)+point(obj.outer_radius_of_slot,0.5*obj.angle_arc_of_tooth)));
            mi_setsegmentprop('<None>',obj.precision_of_stator_iron,obj.automesh,0,obj.group_number_of_stator);
            mi_clearselected;
            
            % draw the slot
            mi_drawarc(point(obj.outer_radius_of_slot,0.5*obj.angle_arc_of_tooth)-[obj.slot_depth, 0],point(obj.inner_radius_of_slot, 0.5*obj.angle_region_of_tooth),0.5*obj.angle_arc_of_auxiliary_slot,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(obj.inner_radius_of_slot,0.5*((0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_auxiliary_slot))+(0.5*obj.angle_region_of_tooth))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawarc(point(obj.outer_radius_of_slot,0.5*(obj.angle_arc_of_tooth)), point(obj.outer_radius_of_slot,0.5*obj.angle_region_of_tooth),0.5*(obj.angle_region_of_tooth-obj.angle_arc_of_tooth),obj.precision_of_stator_iron);
            mi_selectarcsegment(point(obj.outer_radius_of_slot,0.5*((0.5*(obj.angle_arc_of_tooth))+(0.5*obj.angle_region_of_tooth))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator);
            mi_clearselected;
            
            %% Mirror Copy & Rotation Copy
            mi_selectgroup(obj.group_number_of_stator);
            mi_mirror(0,0,obj.outer_radius,0);
            mi_clearselected;
            
            mi_selectgroup(obj.group_number_of_stator);
            mi_copyrotate(obj.base_point,obj.angle_region_of_tooth,( obj.number_of_slots -1));
            mi_clearselected;
            
            %% Draw the Outer radius of VPM
            mi_drawarc(point(obj.outer_radius,0),point(obj.outer_radius,180),180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(obj.outer_radius, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron, '<None>', 0,obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawarc(point(obj.outer_radius,180),point(obj.outer_radius,0),180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(obj.outer_radius, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron, '<None>', 0,obj.group_number_of_stator);
            mi_clearselected;
            %% Draw the Kelvain Transform
            % big circle
            mi_drawarc(point(1.1*obj.outer_radius,0),point(1.1*obj.outer_radius,180),180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(1.1*obj.outer_radius, 0.5*((0)+(180))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron, 'Periodic BC1', 0, obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawarc(point(1.1*obj.outer_radius,180),point(1.1*obj.outer_radius,0),180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(1.1*obj.outer_radius, 0.5*((360)+(180))));
            mi_setarcsegmentprop(obj.precision_of_stator_iron, 'Periodic BC2', 0, obj.group_number_of_stator);
            mi_clearselected;
            
            % small circle
            mi_drawarc(point(0.1*obj.outer_radius,0)+[1.5*obj.outer_radius,0], point(0.1*obj.outer_radius,180)+[1.5*obj.outer_radius,0], 180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(0.1*obj.outer_radius, 0.5*((0)+(180)))+[1.5*obj.outer_radius,0]);
            mi_setarcsegmentprop(obj.precision_of_stator_iron, 'Periodic BC1', 0, obj.group_number_of_stator);
            mi_clearselected;
            
            mi_drawarc(point(0.1*obj.outer_radius,180)+[1.5*obj.outer_radius,0], point(0.1*obj.outer_radius,0)+[1.5*obj.outer_radius,0], 180,obj.precision_of_stator_iron);
            mi_selectarcsegment(point(0.1*obj.outer_radius, 0.5*((360)+(180)))+[1.5*obj.outer_radius,0]);
            mi_setarcsegmentprop(obj.precision_of_stator_iron, 'Periodic BC2', 0, obj.group_number_of_stator);
            mi_clearselected;
            
            mi_addnode(point(1.5*obj.outer_radius,0));
            mi_selectnode(point(1.5*obj.outer_radius,0));
            mi_setnodeprop('A=0', obj.group_number_of_stator);
            mi_clearselected;
            
            %% Set the material of stator
            % set the material of Iron
            stator_Iron_point = point(0.5*(obj.inner_radius+obj.outer_radius),0);
            mi_addblocklabel(stator_Iron_point);
            mi_selectlabel(stator_Iron_point);
            mi_setblockprop('Stator_Iron',obj.automesh,obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator,0);
            mi_clearselected;
            
            % air in the Kelvin Transformation
            air1= point(1.05*obj.outer_radius,90);
            mi_addblocklabel(air1);
            mi_selectlabel(air1);
            mi_setblockprop('Air',obj.automesh,obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator,0);
            mi_clearselected;
            
            air2= point(1.45*obj.outer_radius,0);
            mi_addblocklabel(air2);
            mi_selectlabel(air2);
            mi_setblockprop('Air',obj.automesh,obj.precision_of_stator_iron,'<None>',0,obj.group_number_of_stator,0);
            mi_clearselected;
            
        end
    end
end