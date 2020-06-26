%% Linear subdomain model based on the time harmonics (Version 1.1)
% A general subdomain model which cut off some unnecessary harmonics in the
% motor and save the time;
% Process Oriented Programming
% EEtao
% tangchentao@zju.edu.cn
%
clc; clear;

%% Induce the functions and packages
addpath('.\Analytical_model\Analyical_model_toolbox\');
addpath('.\Analytical_model\Analyical_model_toolbox\time_harmonics\');
addpath('.\Analytical_model\Subdomain_model\time_harmonics\');
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
parameters_of_stator.depth_of_slot = 10e-3;
parameters_of_stator.pole_pairs_of_stator = 5;
parameters_of_stator.number_of_phase = 3;
parameters_of_stator.turns_of_phase = 200;
parameters_of_stator.Jex_amplitude = [0.0e6]; %(A/m^2)

parameters_of_other_part.length_of_stack = 10e-3;
parameters_of_other_part.r_reference = 0.5*(parameters_of_rotor.radius_of_rotor+parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot);
parameters_of_other_part.period_of_machine = gcd(parameters_of_rotor.pole_pairs_of_rotor, parameters_of_stator.number_of_slot);

parameters_of_time_harmonics.current_time_harmonics = []*parameters_of_rotor.pole_pairs_of_rotor;
parameters_of_time_harmonics.PM_time_harmonics = (1:2:5).*parameters_of_rotor.pole_pairs_of_rotor;

parameters_of_space_harmonics.space_harmonics_constant_max = 2;
parameters_of_space_harmonics.lamda_max = 3;

%% Generate the time harmonics and space harmonics collection and determine the unit slot
% time harmonics collection
parameters_of_time_harmonics.time_harmonics_collection = union(parameters_of_time_harmonics.PM_time_harmonics, parameters_of_time_harmonics.current_time_harmonics);
parameters_of_time_harmonics.time_harmonics_collection = parameters_of_time_harmonics.time_harmonics_collection(:)'; % 强制使得时间谐波的合集为行向量
% space harmonics collection
parameters_of_space_harmonics.PM_space_harmonics = PM_space_harmonics(parameters_of_stator, parameters_of_space_harmonics, parameters_of_time_harmonics);
parameters_of_space_harmonics.current_space_harmonics = Current_space_harmonics(parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics);
parameters_of_space_harmonics.space_harmonics_collection = union(parameters_of_space_harmonics.PM_space_harmonics, parameters_of_space_harmonics.current_space_harmonics);
parameters_of_space_harmonics.space_harmonics_collection = parameters_of_space_harmonics.space_harmonics_collection(:)'; % 强制使得空间谐波的合集为行向量
% determine the unit slots 
parameters_of_stator.unit_slots = Unitslots(parameters_of_stator, parameters_of_other_part, parameters_of_time_harmonics);

%% useful coefficients
mu0 = 4*pi*10^-7;
deg = pi/180;

%% Calculate the Winding matrix
parameters_of_windings.number_of_slots = parameters_of_stator.number_of_slot;
parameters_of_windings.pole_pairs_of_stator = parameters_of_stator.pole_pairs_of_stator;
parameters_of_windings.number_of_phase = parameters_of_stator.number_of_phase;
parameters_of_windings.turns_of_phase = parameters_of_stator.turns_of_phase;
parameters_of_windings.pitch_of_coils =  round(parameters_of_windings.number_of_slots/(2*parameters_of_windings.pole_pairs_of_stator));
unit_windingmatrix_of_phaseA= Unit_WindingMatrix_of_PhaseA(parameters_of_windings);
[~, integrated_winding_matrix] = Integrated_WindingMatrix(parameters_of_windings, unit_windingmatrix_of_phaseA);

%% Calculate the Az in the sloted motor
time_harmonics_vertor = parameters_of_time_harmonics.time_harmonics_collection;
space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection;
size_of_time_harmonics_vector = length(time_harmonics_vertor);
size_of_space_harmonics_vector = length(space_harmonics_vector);
size_of_space_harmonics = 4*size_of_space_harmonics_vector+...
    2*parameters_of_stator.unit_slots*(parameters_of_space_harmonics.lamda_max+1);
matrixCNmn = zeros(size_of_time_harmonics_vector, size_of_space_harmonics);

tic;
m = 0;
for time_harmonics = time_harmonics_vertor
    m = m+1;
    matrixPNmn = MatrixPNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, time_harmonics);
    matrixRNmn =  MatrixRNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, time_harmonics, integrated_winding_matrix);
    matrixCNmn_single = linsolve(matrixPNmn, matrixRNmn);
    matrixCNmn(size_of_time_harmonics_vector+m, :) =  matrixCNmn_single';
    matrixCNmn(size_of_time_harmonics_vector-m+1, :) = MatrixCNmn_negative_time_harmonics(parameters_of_space_harmonics, matrixCNmn_single);
end


Omega_mechanical = 600; % (rpm)
omegam = 2*pi*Omega_mechanical/60; % (rad/s)
time_min = 0;
dtime =2*pi/omegam/parameters_of_rotor.pole_pairs_of_rotor/36;
time_max = 2*pi/omegam/parameters_of_rotor.pole_pairs_of_rotor-dtime;
size_of_time = round((time_max-time_min)/dtime+1);
flux_linkage_of_phase = zeros(parameters_of_stator.number_of_phase, size_of_time);

for time = time_min: dtime: time_max
    t = round((time-time_min)/dtime+1);
    thetam = omegam*time;
    matrixAz = MatrixAz_in_the_slots_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, integrated_winding_matrix, matrixCNmn, thetam);
    flux_linkage_of_phase(:, t)= FluxLinkage_basedonAz(parameters_of_stator, parameters_of_other_part, matrixAz, integrated_winding_matrix);
end
toc;

%% Compute the torque
[~, Phi] = FFT(flux_linkage_of_phase(1, :), size_of_time);
phim = Phi(2);
Im = 10;
torque_mean = parameters_of_rotor.pole_pairs_of_rotor*parameters_of_stator.number_of_phase*0.5*phim*Im;

%% Subdomain debug program
% size_of_space_harmonics = 4*(2*parameters_of_space_harmonics.space_harmonics_max)+...
%     2*parameters_of_stator.number_of_slot/parameters_of_other_part.period_of_machine*(parameters_of_space_harmonics.lamda_max+1);
% time_harmonics_vertor = union(parameters_of_time_harmonics.PM_time_harmonics, parameters_of_time_harmonics.current_time_harmonics);
% size_of_time_harmonics = length(time_harmonics_vertor);
% time_harmonics_vertor = union(time_harmonics_vertor, -time_harmonics_vertor);
% matrixCNmn = zeros(size_of_time_harmonics, size_of_space_harmonics);
% m = 0;
% for time_harmonics = time_harmonics_vertor
%     m = m+1;
%     matrixPNmn = MatrixPNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics);
%     matrixRNmn = MatrixRNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, time_harmonics);
%     matrixCNmn(m,:)=  linsolve(matrixPNmn,matrixRNmn);
% end


% r_max = parameters_of_stator.radius_of_stator-0.01e-3;
% r_min = parameters_of_rotor.radius_of_rotor;
% dr = (r_max-r_min)/40;
% size_of_r = round((r_max-r_min)/dr+1);
% phi_min = 0;
% dphi = 2*pi/360;
% phi_max = 2*pi-dphi;
% size_of_phi = round((phi_max-phi_min)/dphi+1);
% for radius = r_min: dr : r_max
%     m = round((radius-r_min)/dr+1);
%     for phi = phi_min: dphi: phi_max
%         n = round((phi-phi_min)/dphi+1);
%         Az(m, n) = Az_subdomain(radius, phi, parameters_of_stator, parameters_of_rotor, parameters_of_other_part, thetam, matrixCN);
%     end
% end
% figure();
% mesh(real(Az));
% 
% thetae_min = 0*pi;
% dthetae = 0.01*pi;
% thetae_max = 2*pi;
% for thetae = thetae_min: dthetae: thetae_max
%     m = round((thetae-thetae_min)/dthetae+1);
%     for j = 1:1:12
%         radius = parameters_of_stator.radius_of_stator+0.5*parameters_of_stator.depth_of_slot;
%         angle = (j-1)*2*pi/parameters_of_stator.number_of_slot+pi/parameters_of_stator.number_of_slot;
%         Az(m,j) = Az_subdomain_time_harmonics(radius, angle, parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, matrixCNmn, thetae);
%     end
% end
% figure();
% plot(real(Az));

% 
% harmonics_max = 50;
% dharmonics = 1;
% size_of_harmonics = round((2*harmonics_max)+1);
% matrixCN = zeros(4, size_of_harmonics);
% for harmonics = -harmonics_max: dharmonics: harmonics_max
%     n = round((harmonics+harmonics_max)/dharmonics+1);
%     if (harmonics~=0)
%         Brn = Bremanence(parameters_of_rotor, harmonics);
%         matrixR = MatrixRn(parameters_of_stator, parameters_of_rotor, harmonics, r_reference);
%         matrixP = MatrixPn(parameters_of_rotor, Brn, harmonics, r_reference, thetam);
%         matrixCN(:, n) =  linsolve(matrixR,matrixP);
%     end
% end
% r_max = parameters_of_stator.radius_of_stator;
% r_min = parameters_of_rotor.radius_of_rotor;
% dr = (r_max-r_min)/40;
% size_of_r = round((r_max-r_min)/dr+1);
% phi_min = 0;
% dphi = 2*pi/360;
% phi_max = 2*pi-dphi;
% size_of_phi = round((phi_max-phi_min)/dphi+1);
% Az = zeros(size_of_r, size_of_phi);
% for radius = r_min: dr : r_max
%     m = round((radius-r_min)/dr+1);
%     for phi = phi_min: dphi: phi_max
%         n = round((phi-phi_min)/dphi+1);
%         Az(m,n) = MagneticVectorPotential(radius, phi, thetam, matrixCN, parameters_of_rotor, parameters_of_stator, harmonics_max, dharmonics, r_reference);
%     end
% end
% 



