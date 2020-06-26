function flux_linkage_matrix_of_phase= FluxLinkage_basedonAz(parameters_of_stator, parameters_of_other_part, matrixAz, integrated_winding_matrix)
%FluxLinkage is used to generator the flux linkage per phase
%% Define some useful parameters
[number_of_phase, number_of_slots_per_phase] = size(integrated_winding_matrix);
number_of_coils_per_phase = round(number_of_slots_per_phase/2);
number_of_slot = round(number_of_phase*number_of_coils_per_phase);
turns = parameters_of_stator.turns_of_phase*2/(number_of_slot*2/number_of_phase); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数
Lstk = parameters_of_other_part.length_of_stack;


%% Calculate the flux linkage per phase
% Calculate the flux linkage of coils
flux_linkage_matrix_of_coils = zeros(number_of_phase, number_of_coils_per_phase);
for m = 1: 1: number_of_phase
    for n = 1: 1: number_of_coils_per_phase
        conductor1 = find(matrixAz(1, :)==abs(integrated_winding_matrix(m, 2*n-1)));
        conductor2 =  find(matrixAz(1, :)==abs(integrated_winding_matrix(m, 2*n)));
        flux_linkage_matrix_of_coils(m, n) = sign(integrated_winding_matrix(m, 2*n-1))*turns*Lstk*matrixAz(2, conductor1)+...
            sign(integrated_winding_matrix(m, 2*n))*turns*Lstk*matrixAz(2, conductor2);
    end
end
% Calculate the flux linkage of phases
% flux_linkage_matrix_of_phase = zeros(number_of_phase,1);
flux_linkage_matrix_of_phase = sum(flux_linkage_matrix_of_coils, 2);

end

