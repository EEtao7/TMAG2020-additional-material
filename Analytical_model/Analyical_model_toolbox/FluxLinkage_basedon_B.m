function flux_linkage_matrix_of_phase= FluxLinkage_basedon_B(parameters_of_stator, parameters_of_other_part, Br_sloted, integrated_winding_matrix)
%FluxLinkage is used to generator the flux linkage per phase
%% Define some useful parameters
r = parameters_of_other_part.radius;
Lstk = parameters_of_other_part.length_of_stack;
turns_of_phase = parameters_of_stator.turns_of_phase;
[number_of_phase, number_of_slots_per_phase] = size(integrated_winding_matrix);
number_of_coils_per_phase = round(number_of_slots_per_phase/2);
number_of_slot = round(number_of_phase*number_of_coils_per_phase);
turns = turns_of_phase*2/(number_of_slot*2/number_of_phase); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数
angle_of_slot = 360/number_of_slot;
dtheta = 360/length(Br_sloted);
deg = pi/180;

%% Calculate the flux linkage per phase
% Calculate the flux linkage of coils
flux_linkage_matrix_of_coils = zeros(number_of_phase, number_of_coils_per_phase);
for m = 1: 1: number_of_phase
    for n = 1: 1: number_of_coils_per_phase
        angle1_of_coil = (abs(integrated_winding_matrix(m, 2*n-1))-1)*angle_of_slot;
        angle2_of_coil = (abs(integrated_winding_matrix(m, 2*n))-1)*angle_of_slot;     
        if(angle1_of_coil<angle2_of_coil)
            theta_min  = angle1_of_coil;
            theta_max = angle2_of_coil;
            for theta = theta_min: dtheta: theta_max-dtheta
                i = round(theta/dtheta+1);
                flux_linkage_matrix_of_coils(m, n) = flux_linkage_matrix_of_coils(m, n)+sign(integrated_winding_matrix(m, 2*n))*turns*Lstk*r*Br_sloted(i)*dtheta*deg;
            end
        elseif(angle1_of_coil>angle2_of_coil)
            theta1 = angle2_of_coil;
            theta2 = angle1_of_coil;
            for theta = 0: dtheta: theta1-dtheta
                i = round(theta/dtheta+1);
                flux_linkage_matrix_of_coils(m, n) = flux_linkage_matrix_of_coils(m, n)+sign(integrated_winding_matrix(m, 2*n))*turns*Lstk*r*Br_sloted(i)*dtheta*deg;
            end
            for theta = theta2: dtheta: 360-dtheta
                i = round(theta/dtheta+1);
                flux_linkage_matrix_of_coils(m, n) = flux_linkage_matrix_of_coils(m, n)+sign(integrated_winding_matrix(m, 2*n))*turns*Lstk*r*Br_sloted(i)*dtheta*deg;
            end
        end
    end
end
% Calculate the flux linkage of phases
flux_linkage_matrix_of_phase = sum(flux_linkage_matrix_of_coils, 2);

end

