function [B_radial, B_tangential] = Magnetic_Field_Noload_Slotless_time_harmonics(parameters_of_rotor, parameters_of_space_harmonics, parameters_of_time_harmonics, matrixBr_PMmn, matrixBt_PMmn, thetam, thetam_shift)
%Magnetic_Field_Noload_Slotless is used to calculate the magnetic field in
% the slotless PMSM without the current (just consider the PM)
%% Define some useful parameters
Pr = parameters_of_rotor.pole_pairs_of_rotor;
time_harmonics_vector = parameters_of_time_harmonics.PM_time_harmonics;
space_harmonics_vector = parameters_of_space_harmonics.PM_space_harmonics;

%% Calculate the magnetic field
% thetam = thetam-thetam_shift; % rotor rotates with mechanical degree shift (time)
B_radial  = 0;
B_tangential = 0;
m = 0;
for time_harmonics = time_harmonics_vector
    m = m+1;
    i = 0;
    for space_harmonics = space_harmonics_vector
        i = i+1;
        n = space_harmonics/Pr;
        % when np ~= 1
        B_radial = B_radial + matrixBr_PMmn(m,i)*cos(n*Pr*thetam-time_harmonics*thetam_shift);
        B_tangential = B_tangential + matrixBt_PMmn(m,i)*sin(n*Pr*thetam-time_harmonics*thetam_shift);
    end
end
B_radial = real(B_radial);
B_tangential = real(B_tangential);

end

