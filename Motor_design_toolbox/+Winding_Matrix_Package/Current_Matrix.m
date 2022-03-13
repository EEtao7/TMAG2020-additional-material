function current_matrix = Current_Matrix(parameters_of_stator, omega, time)
% Current Matrix is used to generate the current matrix for every phase
%% Define some parameters
Np = parameters_of_stator.number_of_phase;

%% Generate the current matrix
angle_offset = 2*pi/Np;
current_matrix = zeros(1, Np);

for m = 1: Np
    current_matrix(1,m) = cos(1*(omega*time-angle_offset*(m-1)));
end

end

