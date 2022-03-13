function lamda = Lamda_function(parameters_of_stator, Lamda0, LamdaNs, thetam)
% This function is used to generate the lamda change with position in the
% air gap.
%% Define some parameters
Ns = parameters_of_stator.number_of_slot;

%% Calculate the relative permeance
lamda = Lamda0-LamdaNs*cos(Ns*thetam);

end

