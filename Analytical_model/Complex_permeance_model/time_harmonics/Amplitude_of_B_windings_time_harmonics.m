function [B_radialn, B_tangentialn] = Amplitude_of_B_windings_time_harmonics(parameters_of_stator, parameters_of_rotor,  parameters_of_other_part, parameters_of_space_harmonics)
%Magnetic_Field_Noload_Slotless is used to calculate the magnetic field in
% the slotless PMSM without the current (just consider the PM)
%% Define some useful Parameters
alpha = parameters_of_stator.ratio_of_tooth;
Rs = parameters_of_stator.radius_of_stator;
Rr = parameters_of_rotor.radius_of_rotor;
Ns = parameters_of_stator.number_of_slot;
% turns_of_phase = parameters_of_stator.turns_of_phase;
r = parameters_of_other_part.radius;
mu0 = pi*4e-7;
b0 = 2*pi/Ns*(1-alpha)*Rs;
% angle_of_slot = 2*pi/Ns;
% [number_of_phase, length] = size(integrated_winding_matrix);
% turns = turns_of_phase*2/(Ns*2/number_of_phase); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数
space_harmonics_vector = parameters_of_space_harmonics.current_space_harmonics;
space_harmonics_vector = space_harmonics_vector(:)';
    
%% Calculate the magnetic field
% current_matrix = zeros(number_of_phase, length);
% for m = 1 : number_of_phase
%     for n = 1 : length
%         current_matrix(m,n) = turns*current_of_every_phase(m)*sign(integrated_winding_matrix(m,n));
%     end
% end

size_of_space_harmonics = length(space_harmonics_vector);
B_radialn  = zeros(1, size_of_space_harmonics);
B_tangentialn = zeros(1, size_of_space_harmonics);

i = 0;
for space_harmonics =space_harmonics_vector
    i = i+1;
    n = space_harmonics;
    B_radialn(i) = mu0/pi/n*(sin(n*b0/2/Rs)/(n*b0/2/Rs))*(n/r*(r/Rs)^n*(1+(Rr/r)^(2*n))/(1-(Rr/r)^(2*n)));
    B_tangentialn(i) = mu0/pi/n*(sin(n*b0/2/Rs)/(n*b0/2/Rs))*(n/r*(r/Rs)^n*(1-(Rr/r)^(2*n))/(1-(Rr/r)^(2*n)));
end
end



