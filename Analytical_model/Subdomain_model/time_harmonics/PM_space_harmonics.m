function n1 = PM_space_harmonics(parameters_of_stator, parameters_of_space_harmonics, parameters_of_time_harmonics)
% This function is used to generate the space harmonics with the PM source
%% Define some useful parameters
c1_max = parameters_of_space_harmonics.space_harmonics_constant_max;
hm = parameters_of_time_harmonics.PM_time_harmonics;
Ns = parameters_of_stator.number_of_slot;
n1 = [];

%% Generate the space harmonics with the PM
for c1 = -c1_max: 1: c1_max
    n_assistant1 = hm+c1*Ns;
    n_assistant2 = -hm+c1*Ns;
    n1 = union(n1, union(n_assistant1, n_assistant2));
end
end

