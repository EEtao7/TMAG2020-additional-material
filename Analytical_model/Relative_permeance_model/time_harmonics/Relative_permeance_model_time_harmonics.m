%% Relative permeance model based on the time harmonics (Version 1.0)
% A general relative permeance model based on the time harmonics
% Process Oriented Programming
% EEtao
% tangchentao@zju.edu.cn
%
clc; clear;

%% Induce the functions and packages
addpath('.\Analytical_model\Analyical_model_toolbox\');
addpath('.\Analytical_model\Relative_permeance_model\time_harmonics\')
addpath('.\Motor_design_toolbox');
addpath('.\FFT_analysis');
import Winding_Matrix_Package.*;

%% Parameters of motor
parameters_of_rotor.B_remanence = 1.2;
parameters_of_rotor.mur = 1.05;
parameters_of_rotor.pole_pairs_of_rotor = 5;
parameters_of_rotor.radius_of_rotor = 53e-3;
parameters_of_rotor.ratio_of_magnet = 1;
parameters_of_rotor.thickness_of_PM = 3e-3;

parameters_of_stator.radius_of_stator = 57e-3;
parameters_of_stator.ratio_of_tooth = 0.8;
parameters_of_stator.number_of_slot = 12;
parameters_of_stator.pole_pairs_of_stator = 5;
parameters_of_stator.number_of_phase = 3;
parameters_of_stator.turns_of_phase = 200;
parameters_of_stator.current_amplitude = []*1;

parameters_of_other_part.length_of_stack = 10e-3;
parameters_of_other_part.period_of_machine = gcd(parameters_of_rotor.pole_pairs_of_rotor, parameters_of_stator.number_of_slot);
parameters_of_other_part.air_gap = parameters_of_stator.radius_of_stator - (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM);
parameters_of_other_part.radius  = parameters_of_stator.radius_of_stator - 0.5*parameters_of_other_part.air_gap;

% time harmonics collection
parameters_of_time_harmonics.current_time_harmonics = []*parameters_of_rotor.pole_pairs_of_rotor;
parameters_of_space_harmonics.space_harmonics_constant_max = 2;
parameters_of_time_harmonics.PM_time_harmonics = (1:2:5).*parameters_of_rotor.pole_pairs_of_rotor;

% space harmonics collection
parameters_of_space_harmonics.PM_space_harmonics = parameters_of_time_harmonics.PM_time_harmonics;
parameters_of_space_harmonics.current_space_harmonics = Current_space_harmonics(parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics);

%% Other parameters
deg = pi/180;
thetam_min = 0;
dthetam = 360/360;
thetam_max = 360-dthetam;
size_of_thetam = round((thetam_max-thetam_min)/dthetam+1);

tic;
%% Calculate the relative permeance function lamda in the motor
lamda = zeros(1, size_of_thetam);
[Lamda0, LamdaNs] = Relative_Permeance_conformal_transformation (parameters_of_stator, parameters_of_rotor, parameters_of_other_part);
for thetam = thetam_min : dthetam : thetam_max
    m = round((thetam-thetam_min)/dthetam+1);
    lamda(m) = Lamda_function(parameters_of_stator, Lamda0, LamdaNs, thetam*deg);
end

%% Calculate the Winding matrix
parameters_of_windings.number_of_slots = parameters_of_stator.number_of_slot;
parameters_of_windings.pole_pairs_of_stator = parameters_of_stator.pole_pairs_of_stator;
parameters_of_windings.number_of_phase = parameters_of_stator.number_of_phase;
parameters_of_windings.turns_of_phase = parameters_of_stator.turns_of_phase;
parameters_of_windings.pitch_of_coils =  round(parameters_of_windings.number_of_slots/(2*parameters_of_windings.pole_pairs_of_stator));
unit_windingmatrix_of_phaseA= Unit_WindingMatrix_of_PhaseA(parameters_of_windings);
[~, integrated_winding_matrix] = Integrated_WindingMatrix(parameters_of_windings, unit_windingmatrix_of_phaseA);

%% Compute the PM coefficient matrix
time_harmonics_vector = parameters_of_time_harmonics.PM_time_harmonics;
space_harmonics_vector = parameters_of_space_harmonics.PM_space_harmonics;
size_of_time_harmonics_vector = length(time_harmonics_vector);
size_of_space_harmonics_vector = length(space_harmonics_vector);
matrixBr_PMmn = zeros(size_of_time_harmonics_vector, size_of_space_harmonics_vector);
matrixBt_PMmn = zeros(size_of_time_harmonics_vector, size_of_space_harmonics_vector);
m = 0;
for time_harmonics = time_harmonics_vector
    m = m+1;
    [matrixBr_PMn, matrixBt_PMn] = Amplitude_of_B_PM_time_harmonics(parameters_of_stator, parameters_of_rotor,  parameters_of_other_part, parameters_of_space_harmonics, time_harmonics);
     matrixBr_PMmn(m, :) = matrixBr_PMn;
    matrixBt_PMmn(m, :) = matrixBt_PMn;
end

%% Rotate the Magnetic with the time
omegam = 2*pi*50; % mechanical speed;
time_min = 0;
dtime =2*pi/parameters_of_rotor.pole_pairs_of_rotor/omegam/36;
time_max = 2*pi/parameters_of_rotor.pole_pairs_of_rotor/omegam-dtime;
size_of_time = round((time_max-time_min)/dtime+1);
B_radial_slotless_PM = zeros(1, size_of_thetam);
B_tangential_slotless_PM = zeros(1, size_of_thetam);
B_radial_slotless_windings = zeros(1, size_of_thetam);
B_tangential_slotless_windings = zeros(1, size_of_thetam);
B_radial_slotless = zeros(1, size_of_thetam);
B_tangential_slotless = zeros(1, size_of_thetam);
B_radial_sloted = zeros(1, size_of_thetam);
B_tangential_sloted = zeros(1, size_of_thetam);
flux_linkage_of_phase = zeros(parameters_of_stator.number_of_phase, size_of_time);

for time = time_min: dtime: time_max
    t = round((time-time_min)/dtime+1);
    thetam_shift = omegam*time; % mechanical degree;
    current_of_every_phase = Current_of_every_phase(parameters_of_stator, parameters_of_rotor, parameters_of_time_harmonics, thetam_shift);
    for thetam = thetam_min : dthetam : thetam_max
        m = round((thetam-thetam_min)/dthetam+1);
        % Compute the magnetic field generated by PM
        [B_radial_slotless_PM(m), ~] = Magnetic_Field_Noload_Slotless_time_harmonics(parameters_of_rotor,...
            parameters_of_space_harmonics, parameters_of_time_harmonics, matrixBr_PMmn, matrixBt_PMmn, thetam*deg, thetam_shift);
        % Compute the magnetic field generated by windings
        [B_radial_slotless_windings(m), ~] = Magnetic_Field_of_Current_time_harmonics( parameters_of_stator,...
            parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, current_of_every_phase, integrated_winding_matrix, thetam*deg);
        % Compute the total magnetic field by superposition
        B_radial_slotless(m) = B_radial_slotless_PM(m)+B_radial_slotless_windings(m);
        % Calculate the Magnetic field in the Sloted Motor
         B_radial_sloted(m) = B_radial_slotless(m)*lamda(m);
    end
  % Calculate the flux linkage of the windings
    flux_linkage_of_phase(:, t)= FluxLinkage_basedon_B(parameters_of_stator, parameters_of_other_part, B_radial_sloted, integrated_winding_matrix);
end
toc;

%% Calculate the electrical torque
[~, Phi] = FFT(flux_linkage_of_phase(1, :), size_of_time);
phim = Phi(2);
Im = 10;
torque_mean = parameters_of_rotor.pole_pairs_of_rotor*parameters_of_stator.number_of_phase*0.5*phim*Im;

