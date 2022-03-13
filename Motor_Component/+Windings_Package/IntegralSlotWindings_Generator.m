classdef IntegralSlotWindings_Generator< Windings_Package.WindingMatrix_Generator
    % define the GeneralFractionalSlotsWindings class
    properties
        Winding_Matrix_PhaseA
        turns_of_windings
    end
    
    %% define the Dependent properties
    properties(Dependent)
        machine_period
        pitch_of_windings
        slots_per_phase_per_pole
        number_of_slots_per_unit
    end
    
    methods
        %% Dependent properties caculation function
        function machine_period = get.machine_period(obj)
            machine_period = gcd(obj.number_of_slots, obj.pole_pairs_of_windings);
        end
        
        function slots_per_phase_per_pole = get.slots_per_phase_per_pole(obj)
            %             slots_per_phase_per_pole = obj.number_of_slots*obj.number_of_layers/obj.number_of_phase/(2*obj.pole_pairs_of_windings);
            slots_per_phase_per_pole = obj.number_of_slots/obj.number_of_phase/(2*obj.pole_pairs_of_windings);
        end
        
        function number_of_slots_per_unit = get.number_of_slots_per_unit(obj)
            number_of_slots_per_unit = obj.number_of_slots * obj.number_of_layers / obj.number_of_phase / obj.machine_period;
        end
        
        function pitch_of_windings = get.pitch_of_windings(obj)
            if(obj.number_of_layers == 1)
                pitch_of_windings = obj.number_of_slots/(2*obj.pole_pairs_of_windings);
            elseif(obj.number_of_layers == 2)
                pitch_of_windings = obj.number_of_slots/(2*obj.pole_pairs_of_windings)-1;
            end
        end
        
        %% Various function of femm
        function define_name(obj)
            %define the name of windings
            obj.name = strcat('ISW_', num2str(obj.number_of_slots),'-', num2str(obj.number_of_layers),...
                '-', num2str(obj.pole_pairs_of_windings),'-', num2str(obj.slots_per_phase_per_pole));
        end
        
        function generate_the_windingmatrix_phaseA(obj)
            % generate the winding matrix of phaseA
            % initialize the winding matrix phaseA
            obj.Winding_Matrix_PhaseA = zeros(3, obj.number_of_slots_per_unit);
            
            % add the element in the winding matrix phaseA
            % for the one layer windings
           if(obj.number_of_layers == 1)
                if(obj.pitch_of_windings == obj.number_of_slots/(2*obj.pole_pairs_of_windings))
                    for m = 1:1:round(obj.slots_per_phase_per_pole)
                        obj.Winding_Matrix_PhaseA(1,m)=m;
                        obj.Winding_Matrix_PhaseA(2,m)=1;
                        obj.Winding_Matrix_PhaseA(3,m)=obj.turns_of_windings;
                        obj.Winding_Matrix_PhaseA(1,m+round(obj.slots_per_phase_per_pole))=m+obj.pitch_of_windings;
                        obj.Winding_Matrix_PhaseA(2,m+round(obj.slots_per_phase_per_pole))=1;
                        obj.Winding_Matrix_PhaseA(3,m+round(obj.slots_per_phase_per_pole))=-obj.turns_of_windings;
                    end
                else
                    print('pitch of windings have errors');
                end
                % for the two layers windings
            elseif(obj.number_of_layers == 2)
                for m = 1:1:round(obj.slots_per_phase_per_pole)
                    obj.Winding_Matrix_PhaseA(1,m)=m;
                    obj.Winding_Matrix_PhaseA(2,m)=2; % 上层绕组为2
                    obj.Winding_Matrix_PhaseA(3,m)=obj.turns_of_windings;
                    obj.Winding_Matrix_PhaseA(1,round(m+obj.slots_per_phase_per_pole))=m+obj.pitch_of_windings+1;
                    obj.Winding_Matrix_PhaseA(2,round(m+obj.slots_per_phase_per_pole))=2; % 上层绕组为2
                    obj.Winding_Matrix_PhaseA(3,round(m+obj.slots_per_phase_per_pole))=-obj.turns_of_windings;
                    
                    obj.Winding_Matrix_PhaseA(1,round(m+2*obj.slots_per_phase_per_pole))=m-1;
                    obj.Winding_Matrix_PhaseA(2,round(m+2*obj.slots_per_phase_per_pole))=1; % 下层绕组为1
                    obj.Winding_Matrix_PhaseA(3,round(m+2*obj.slots_per_phase_per_pole))=obj.turns_of_windings;
                    obj.Winding_Matrix_PhaseA(1,round(m+3*obj.slots_per_phase_per_pole))=m+obj.pitch_of_windings;
                    obj.Winding_Matrix_PhaseA(2,round(m+3*obj.slots_per_phase_per_pole))=1; % 下层绕组为1
                    obj.Winding_Matrix_PhaseA(3,round(m+3*obj.slots_per_phase_per_pole))=-obj.turns_of_windings;
                end
                %check the winding matrix
                for m = 1:1:obj.number_of_slots_per_unit
                    if(obj.Winding_Matrix_PhaseA(1,m) <= 0)
                        obj.Winding_Matrix_PhaseA(1,m) = obj.Winding_Matrix_PhaseA(1,m)+obj.number_of_slots;
                    elseif(obj.Winding_Matrix_PhaseA(1,m) > obj.number_of_slots_per_unit)
                        obj.Winding_Matrix_PhaseA(1,m) = obj.Winding_Matrix_PhaseA(1,m)-obj.number_of_slots;
                    end
                end
            else
                print('Generator can not generate the winding matrix ');
            end
            
        end
    end
end