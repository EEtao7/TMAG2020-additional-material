function [Lamda0, LamdaNs] = Relative_Permeance_conformal_transformation (parameters_of_stator, parameters_of_rotor, parameters_of_other_part)
% This function is used to compute the relative permeance coefficient 
%% Define some parameters
length_of_air_gap = parameters_of_other_part.air_gap;
thickness_of_PM = parameters_of_rotor.thickness_of_PM;
mur = parameters_of_rotor.mur;
Ns = parameters_of_stator.number_of_slot;
Rs = parameters_of_stator.radius_of_stator;
alpha = parameters_of_stator.ratio_of_tooth;
pitch_of_tooth = Rs*2*pi/Ns;
slot = Rs*2*pi/Ns*(1-alpha);

equivalent_airgap = length_of_air_gap+thickness_of_PM/mur;
x = slot/pitch_of_tooth;
y = pitch_of_tooth/equivalent_airgap;
beta = 0.5-0.5/sqrt(1+(0.5*x*y)^2);

%% Calculate the relative permeance
Lamda0 = 1-1.6*beta*x;
LamdaNs = 4/pi*beta*(0.5+x^2/(0.78125-2*x^2))*sin(1.6*pi*x);

end

