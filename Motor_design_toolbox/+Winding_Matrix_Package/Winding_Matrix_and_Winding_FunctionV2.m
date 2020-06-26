%% Winding matrix and Winding function Version (2.0)
%
% Processed Programming
% EEtao
% tangchentao@zju.edu.cn
%
clc; clear;

%% Induce the functions and packages
addpath('.\Motor_design_toolbox');
addpath('.\FFT_analysis');
import Winding_Matrix_Package.*;

%% Define some parameters
parameters_of_windings.number_of_slots = 54;
parameters_of_windings.pole_pairs_of_stator = 3;
parameters_of_windings.number_of_phase = 3;
parameters_of_windings.turns_of_phase = 200;
parameters_of_windings.pitch_of_coils =  round(parameters_of_windings.number_of_slots/(2*parameters_of_windings.pole_pairs_of_stator));

thetam.thetam_min = 0;
thetam.dthetam = 0.5;
thetam.thetam_max = 360-thetam.dthetam;
thetam.size_of_thetam = round((thetam.thetam_max - thetam.thetam_min)/ thetam.dthetam+1);

%% Calculate the Winding matrix
% Calculate the Winding matrix of phaseA in the unit machine
unit_windingmatrix_of_phaseA= Unit_WindingMatrix_of_PhaseA(parameters_of_windings);
% Calculate the Winding matrix of PhaseA and Winding matrix for the whole machine
[integrated_windingmatrix_of_phaseA, integrated_winding_matrix] = Integrated_WindingMatrix(parameters_of_windings, unit_windingmatrix_of_phaseA);

%% Calculate the winding function
% winding function of phaseA
winding_function_PhaseA = Winding_function_of_phase(parameters_of_windings, thetam, integrated_windingmatrix_of_phaseA);
% winding function for every phase
winding_function = zeros(parameters_of_windings.number_of_phase, thetam.size_of_thetam);
for m = 1: 1: parameters_of_windings.number_of_phase
    winding_function(m, :) = Winding_function_of_phase(parameters_of_windings, thetam, integrated_winding_matrix(m, :));
end

%% Calculate the winding factor based on the winding function
[pole_pairs_of_MMF, windingfactormatrix] = WindingFactorMatrix(parameters_of_windings, thetam, winding_function_PhaseA);

%% Calculate the MMF of whole windings
% Generate the current matrix
omega = 2*pi*50;
time = 0;
current_matrix = Current_Matrix(parameters_of_windings, omega, time);

% Compute the MMF
MMF = current_matrix*winding_function;
[f_MMF,P_MMF]=FFT(MMF, thetam.size_of_thetam);

% Draw the figure of MMF and Winding factor
figure();
subplot(3,1,1);
plot(thetam.thetam_min: thetam.dthetam: thetam.thetam_max, MMF);
xlabel('Mechanical Angle, Degrees');
ylabel('Amplitude, A*Turns');
title('MMF of the windings');
grid on;
subplot(3,1,2);
bar(f_MMF, P_MMF, 1.0);
xlabel('Harmonics Order');
ylabel('Amplitude, A*Turns');
title('MMF of the windings (FFT)');
grid on;
subplot(3,1,3);
bar(pole_pairs_of_MMF, windingfactormatrix, 0.4);
xlabel('Harmonics Order');
ylabel('Amplitude');
title('Winding factor of the windings');
grid on;

