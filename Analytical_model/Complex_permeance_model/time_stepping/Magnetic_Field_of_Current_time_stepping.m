function[B_radial, B_tangential] = Magnetic_Field_of_Current_time_stepping( parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, current_of_every_phase, integrated_winding_matrix, thetam)
% Magnetic_Field_of_Current is used to calculate the magnetic field in the
% slotless PMSM without the PM (just considet the current)
%% Define some useful Parameters
alpha = parameters_of_stator.ratio_of_tooth;
Rs = parameters_of_stator.radius_of_stator;
Rr = parameters_of_rotor.radius_of_rotor;
Ns = parameters_of_stator.number_of_slot;
turns_of_phase = parameters_of_stator.turns_of_phase;
r = parameters_of_other_part.radius;
mu0 = pi*4e-7;
b0 = 2*pi/Ns*(1-alpha)*Rs;
angle_of_slot = 2*pi/Ns;
[number_of_phase, length] = size(integrated_winding_matrix);
turns = turns_of_phase*2/(Ns*2/number_of_phase); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数
harmonics_max = parameters_of_space_harmonics.windings_space_harmonics_max;

%% Calculate the magnetic field
current_matrix = zeros(number_of_phase, length);
for m = 1 : number_of_phase
    for n = 1 : length
        current_matrix(m,n) = turns*current_of_every_phase(m)*sign(integrated_winding_matrix(m,n));
    end
end

B_radial  = 0;
B_tangential = 0;
for k = 1: 1 : round(number_of_phase*length)
    for n =1: 1: harmonics_max
        B_radial = B_radial + mu0/pi*current_matrix(k)/n*(sin(n*b0/2/Rs)/(n*b0/2/Rs))*(n/r*(r/Rs)^n*(1+(Rr/r)^(2*n))/(1-(Rr/r)^(2*n)))*sin(n*(thetam-(abs(integrated_winding_matrix(k))-1)*angle_of_slot));
        B_tangential = B_tangential + mu0/pi*current_matrix(k)/n*(sin(n*b0/2/Rs)/(n*b0/2/Rs))*(n/r*(r/Rs)^n*(1-(Rr/r)^(2*n))/(1-(Rr/r)^(2*n)))*cos(n*(thetam-(abs(integrated_winding_matrix(k))-1)*angle_of_slot));
    end
end
end

