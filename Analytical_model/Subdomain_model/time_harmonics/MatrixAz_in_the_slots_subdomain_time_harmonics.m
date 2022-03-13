function matrixAz = MatrixAz_in_the_slots_subdomain_time_harmonics(parameters_of_stator,  parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, integrated_winding_matrix, matrixCNmn, thetam)
% MatrixAz is used to generate the Az in the slots based on the subdomain when motor rotating
%% Defime some useful parameters
Ns = parameters_of_stator.number_of_slot;
radius = parameters_of_stator.radius_of_stator+0.5*parameters_of_stator.depth_of_slot;

%% Calculate the magnetic vector potential based on the subdomain
% Compute the matrixAz
matrixAz = zeros(2, Ns);
matrixAz(1,:) = 1: 1: Ns;
for j = 1: 1: Ns  
    angle = (j-1)*2*pi/Ns+pi/Ns;
    Az = Az_subdomain_time_harmonics(radius, angle, parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, integrated_winding_matrix, matrixCNmn, thetam);
    matrixAz(2,j) = real(Az);
end

end

