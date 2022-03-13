function [pole_pairs_of_MMF, windingfactormatrix] = WindingFactorMatrix(parameters_of_stator, thetam, winding_function_of_phase)
%WindingFactorMatrix is used to generate the winding factor matrix based on
%the winding function of one phase
%% Define some parameters 
Ns = parameters_of_stator.number_of_slots;
Ps = parameters_of_stator.pole_pairs_of_stator;
t = gcd(Ns, Ps);
turns_of_phase = parameters_of_stator.turns_of_phase;
size_of_FFT = thetam.size_of_thetam;

%% Calculate the winding factor 
% the low frequency part has some errors of the FFT. For this reason, we need to use the FFT V2 based on the base frequency
[frequency, p_assistant]=FFT_V2(winding_function_of_phase, size_of_FFT, t, 1, Ns); 
pole_pairs_of_MMF = frequency(2: end);
P = p_assistant(2 : end);
windingfactormatrix = zeros(1, Ns/t); % the size of winding factor matrix is the number of slot
for m = 1: floor(0.5*Ns/t)
    windingfactormatrix(1, m) = P(m)/(2/(pi*(m*t))*turns_of_phase);
    windingfactormatrix(1, end-m) = windingfactormatrix(1, m);
end
end

