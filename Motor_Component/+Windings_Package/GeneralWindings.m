classdef GeneralWindings< Windings_Package.Windings
    % define the GeneralFractionalSlotsWindings class
    % This function has the limitation of number of phase equal to 3
    properties
        Winding_Matrix_PhaseA 
        series_flag
        turns_of_phase
    end
    
    %% define the Dependent properties
    properties(Dependent)
        number_of_coils
        machine_period
        number_of_slots_per_unit
    end
    
    methods
        %         %% Constructor function
        %         function obj = ConcentratedWindings(name, precision_of_windings, group_number_of_windings, number_of_phase,...
        %                series_flag, turns, current_of_windings, stator)
        %             obj = obj@Windings(name, precision_of_windings, group_number_of_windings, number_of_phase);
        %             obj.series =series_flag;
        %             obj.turns = turns;
        %             obj.current_of_windings = current_of_windings;
        %             obj.stator = stator;
        %         end
        
        %% Dependent properties caculation function
        function number_of_coils = get.number_of_coils(obj)
            number_of_coils = obj.stator.number_of_slots*obj.number_of_layers/2;
        end
        
        function machine_period = get.machine_period(obj)
            %             if (obj.number_of_layers==1)
            %                 machine_period = gcd(obj.stator.number_of_slots/2, obj.stator.pole_of_pairs_stator);
            %             else
            %                 machine_period = gcd(obj.stator.number_of_slots, obj.stator.pole_of_pairs_stator);
            %             end
            machine_period = gcd(obj.stator.number_of_slots, obj.stator.pole_pairs_of_stator);
        end
        %         % 半齿作为一种特殊的全齿结构可以使用这种统一的判定的方式，但绕组矩阵也需要做出一定的调整
        
        function number_of_slots_per_unit = get.number_of_slots_per_unit(obj)
            number_of_slots_per_unit = obj.stator.number_of_slots*obj.number_of_layers/obj.number_of_phase/obj.machine_period;
        end
        
        %% Various function of femm
        function define_material(obj)
            % Properties of Material
            % define the properties of Coils
            mi_addmaterial('Coil', 1, 1, 0, 0, 58*0.65, 0, 0, 1, 0, 0, 0, 0, 0);
            % define the rank of windings
            mi_addcircprop('A', obj.current_of_windings.ia, obj.series_flag);
            mi_addcircprop('B', obj.current_of_windings.ib, obj.series_flag);
            mi_addcircprop('C', obj.current_of_windings.ic, obj.series_flag);
        end
        
        function draw_single_element(obj)        
            %% divide the slot into different layers
            %define some useful parameters    
            alpha1 = obj.stator.angle_region_of_tooth;
            alpha2 = obj.stator.angle_arc_of_tooth;
            r1 = obj.stator.inner_radius_of_slot;
            r2 = obj.stator.outer_radius_of_slot;
            
            % define some useful radius and angles in different regions
            switch(obj.number_of_layers)
                case 1
                    % The number of layers equal to 1
                    radius_region1 = 0.5*(r1+r2);
                    angle_region1 = 0.5*alpha1;
                case 2
                    % The number of layers equal to 2
                    radius_region1 = 0.5*(r1+r2);
                    angle_region1 = 0.25*alpha1+0.25*alpha2;
                    radius_region2 = 0.5*(r1+r2);
                    angle_region2= 0.75*alpha1-0.25*alpha2;
                    mi_drawline(point(obj.stator.outer_radius_of_slot, 0.5*obj.stator.angle_region_of_tooth), point(obj.stator.inner_radius_of_slot,0.5*obj.stator.angle_region_of_tooth));
                    mi_selectsegment(0.5*(point(obj.stator.outer_radius_of_slot, 0.5*obj.stator.angle_region_of_tooth)+point(obj.stator.inner_radius_of_slot,0.5*obj.stator.angle_region_of_tooth)));
                    mi_setsegmentprop('<None>', obj.stator.precision_of_stator_iron, obj.stator.automesh, 0, obj.stator.group_number_of_stator);
                    mi_copyrotate(obj.stator.base_point, obj.stator.angle_region_of_tooth, ( obj.stator.number_of_slots -1));
                    mi_clearselected;
                case 4
                    % The number of layers equal to 4;
                    radius_region1 = 0.25*r1+0.75*r2;
                    angle_region1 = 0.25*alpha1+0.25*alpha2;
                    radius_region2 = 0.75*r1+0.25*r2;
                    angle_region2 = 0.25*alpha1+0.25*alpha2;
                    radius_region3 = 0.25*r1+0.75*r2;
                    angle_region3 = 0.75*alpha1-0.25*alpha2;
                    radius_region4 = 0.75*r1+0.25*r2;
                    angle_region4 = 0.75*alpha1-0.25*alpha2;
                    mi_drawline(point(obj.stator.outer_radius_of_slot, 0.5*obj.stator.angle_region_of_tooth), point(obj.stator.inner_radius_of_slot,0.5*obj.stator.angle_region_of_tooth));
                    mi_selectsegment(0.5*(point(obj.stator.outer_radius_of_slot, 0.5*obj.stator.angle_region_of_tooth)+point(obj.stator.inner_radius_of_slot,0.5*obj.stator.angle_region_of_tooth)));
                    mi_setsegmentprop('<None>', obj.stator.precision_of_stator_iron, obj.stator.automesh, 0, obj.group_number_of_winding);
                    mi_clearselected;
                    mi_drawline(point(obj.stator.inner_radius_of_slot, 0.5*obj.stator.angle_arc_of_tooth)+[0.5*obj.stator.slot_depth, 0], point(0.5*(obj.stator.outer_radius_of_slot+obj.stator.inner_radius_of_slot),0.5*obj.stator.angle_region_of_tooth));
                    mi_selectsegment(0.5*(point(obj.stator.inner_radius_of_slot, 0.5*obj.stator.angle_arc_of_tooth)+[0.5*obj.stator.slot_depth, 0]+point(0.5*(obj.stator.outer_radius_of_slot+obj.stator.inner_radius_of_slot),0.5*obj.stator.angle_region_of_tooth)));
                    mi_setsegmentprop('<None>', obj.stator.precision_of_stator_iron, obj.stator.automesh, 0, obj.group_number_of_winding);
                    x= point(10*obj.stator.outer_radious, 0.5*obj.stator.angle_region_of_tooth);
                    mi_mirror(0, 0, x(1),x(2));
                    mi_selectgroup(obj.group_number_of_winding);
                    mi_copyrotate(obj.stator.base_point, obj.stator.angle_region_of_tooth, ( obj.stator.number_of_slots -1));
                    mi_clearselected;
                    
            end
            
            %% Draw phase A windings
            %check the slots and poles combination
            if(mod(obj.stator.number_of_slots/obj.machine_period,obj.number_of_phase)==0)
                disp('The Windings satify the slots/poles combination');
            end
            % check the winding matrix of PhaseA
            [~,N] = size(obj.Winding_Matrix_PhaseA);
            if( N ~= obj.number_of_slots_per_unit)
                disp('Winding Matrix PhaseA have errors');
            end
            % define some useful parameters
            angle_of_slot = 360/(obj.stator.number_of_slots);
            
            % choose the layers
            switch(obj.number_of_layers)
                case 1 %number of layer = 1
                    % draw the phaseA coils
                    for n = 1:N
                        coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1);
                        mi_addblocklabel(coil_point);
                        mi_selectlabel(coil_point);
                        mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                        mi_clearselected;
                    end
                    mi_selectgroup(obj.group_number_of_windingA);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                    mi_clearselected;
                    
                case 2 %number of layer = 2
                    for n = 1:N
                        % draw the phaseA coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingA);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
                    
                case 4 %number of layer = 4
                    for n = 1:N
                        % draw the phaseA coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 3
                                coil_point = point(radius_region3, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region3);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 4
                                coil_point = point(radius_region4, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region4);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'A', 0, obj.group_number_of_windingA, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingA);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
            end
            
            % Draw phase B windings
            % check the winding matrix of PhaseA
            [~,N] = size(obj.Winding_Matrix_PhaseA);
            if( N ~= obj.number_of_slots_per_unit)
                disp('Winding Matrix PhaseA have errors');
            end
            % define some useful parameters
            angle_of_windingB = 120/obj.machine_period;% 由于绕组分布理论基于星形图理论，所有可以发现电机周期这个概念要比绕组极对数更加具有普遍性
            % choose the layers
            switch(obj.number_of_layers)
                case 1 %number of layer = 1
                    % draw the phaseB coils
                    for n = 1:N
                        coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingB);
                        mi_addblocklabel(coil_point);
                        mi_selectlabel(coil_point);
                        mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                        mi_clearselected;
                    end
                    mi_selectgroup(obj.group_number_of_windingB);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
                    
                case 2 %number of layer = 2
                    for n = 1:N
                        % draw the phaseB coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingB);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
                    
                case 4 %number of layer = 4
                    for n = 1:N
                        % draw the phaseB coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 3
                                coil_point = point(radius_region3, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region3+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 4
                                coil_point = point(radius_region4, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region4+angle_of_windingB);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'B', 0, obj.group_number_of_windingB, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingB);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
            end
            
            % Draw phase C windings
            % check the winding matrix of PhaseA
            [~,N] = size(obj.Winding_Matrix_PhaseA);
            if( N ~= obj.number_of_slots_per_unit)
                disp('Winding Matrix PhaseA have errors');
            end
            % define some useful parameters
            angle_of_windingC = -120/obj.machine_period;
            % choose the layers
            switch(obj.number_of_layers)
                case 1 %number of layer = 1
                    % draw the phaseC coils
                    for n = 1:N
                        coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingC);
                        mi_addblocklabel(coil_point);
                        mi_selectlabel(coil_point);
                        mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                        mi_clearselected;
                    end
                    mi_selectgroup(obj.group_number_of_windingC);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
                    
                case 2 %number of layer = 2
                    for n = 1:N
                        % draw the phaseC coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingC);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
                    
                case 4 %number of layer = 4
                    for n = 1:N
                        % draw the phaseC coils
                        switch (obj.Winding_Matrix_PhaseA(2,n))
                            %choose the regions
                            case 1
                                coil_point = point(radius_region1, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region1+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 2
                                coil_point = point(radius_region2, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region2+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 3
                                coil_point = point(radius_region3, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region3+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                            case 4
                                coil_point = point(radius_region4, (obj.Winding_Matrix_PhaseA(1,n)-1)*angle_of_slot+angle_region4+angle_of_windingC);
                                mi_addblocklabel(coil_point);
                                mi_selectlabel(coil_point);
                                mi_setblockprop('Coil', obj.automesh, obj.precision_of_windings,'C', 0, obj.group_number_of_windingC, obj.Winding_Matrix_PhaseA(3,n));
                                mi_clearselected;
                        end
                    end
                    mi_selectgroup(obj.group_number_of_windingC);
                    mi_copyrotate(obj.base_point,360/obj.machine_period, obj.machine_period-1);
                    mi_clearselected;
            end
            
        end
    end
end