%% Flux linkage Test Based on Linear Subdomain Model  (Vernier PMSM)
% This program is designed to test the flux linkage in the Vernier PMSM
% based on the linear subdomain model
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

%% Initial Parameters of motor
% parameters of rotor
parameters_of_rotor.B_remanence = 1.2;
parameters_of_rotor.mur = 1.05;
parameters_of_rotor.pole_pairs_of_rotor = 5;
parameters_of_rotor.radius_of_rotor = 53e-3;
parameters_of_rotor.ratio_of_magnet = 1;
parameters_of_rotor.thickness_of_PM = 3e-3;

% parameters of stator
parameters_of_stator.radius_of_stator = 57e-3;
parameters_of_stator.ratio_of_tooth = 0.9;
parameters_of_stator.number_of_slot = 12;
parameters_of_stator.depth_of_slot = 20e-3;
parameters_of_stator.pole_pairs_of_stator = 5;
parameters_of_stator.number_of_phase = 3;
parameters_of_stator.turns_of_phase = 200;
parameters_of_stator.Jex_amplitude = 0; %(A/m^2)

% parameters of other part
parameters_of_other_part.length_of_stack = 10e-3;

% uniformize the geometry parameters
parameters_of_other_part.unit_coefficient = 0.5*(parameters_of_rotor.radius_of_rotor+parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot); 
parameters_of_other_part.r_reference1 = 0.5*(parameters_of_rotor.radius_of_rotor+parameters_of_stator.radius_of_stator);
parameters_of_other_part.r_reference2 = 0.5*(parameters_of_stator.radius_of_stator+(parameters_of_stator.radius_of_stator++parameters_of_stator.depth_of_slot));
parameters_of_other_part.period_of_machine = gcd(parameters_of_rotor.pole_pairs_of_rotor, parameters_of_stator.number_of_slot);

% parameters of time harmonics
parameters_of_time_harmonics.current_time_harmonics = []*parameters_of_rotor.pole_pairs_of_rotor;
parameters_of_time_harmonics.PM_time_harmonics = [1].*parameters_of_rotor.pole_pairs_of_rotor;

% parameters of space harmonics
parameters_of_space_harmonics.space_harmonics_constant_max = 2;
parameters_of_space_harmonics.lamda_max = 3;

% parameters of performance
parameters_of_performance.flux_linkage = 0;

%% Calculate the flux linkage in the different motors
% Range of  parameters which needs to be adjusted
% number of slots
Ns_max = 8;
Ns_min = 3;
dNs = 1;
Ns_range = parameters_of_stator.number_of_phase*(Ns_min: dNs: Ns_max);

% pole-pairs of the parameters
Pr_max = 23;
Pr_min = 8;
dPr = 1;
Pr_range = Pr_min: dPr: Pr_max;

% ratio of tooth
Ratio_of_tooth_max = 0.9;
Ratio_of_tooth_min = 0.3;
dRatio_of_tooth = 0.1;
Ratio_of_tooth_range = Ratio_of_tooth_min: dRatio_of_tooth: Ratio_of_tooth_max;

% Calculate and save the data in the Subdomain model
tic
n = 0;
for Ns = Ns_range
    parameters_of_stator.number_of_slot = Ns;
    for Pr = Pr_range
        parameters_of_rotor.pole_pairs_of_rotor = Pr;
        Ps = Ns-Pr; % for the Vernier PMSM
        parameters_of_stator.pole_pairs_of_stator = Ps; 
        for Ratio_of_tooth = Ratio_of_tooth_range
            parameters_of_stator.ratio_of_tooth = Ratio_of_tooth;
             %% adjust some assistant parameters with the key parameters
            parameters_of_other_part.unit_coefficient = 0.5*(parameters_of_rotor.radius_of_rotor+parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot);
            parameters_of_other_part.r_reference1 = 0.5*(parameters_of_rotor.radius_of_rotor+parameters_of_stator.radius_of_stator);
            parameters_of_other_part.r_reference2 = 0.5*(parameters_of_stator.radius_of_stator+(parameters_of_stator.radius_of_stator++parameters_of_stator.depth_of_slot));
            parameters_of_other_part.period_of_machine = gcd(parameters_of_rotor.pole_pairs_of_rotor, parameters_of_stator.number_of_slot);
            parameters_of_time_harmonics.current_time_harmonics = []*parameters_of_rotor.pole_pairs_of_rotor;
            parameters_of_time_harmonics.PM_time_harmonics = (1:2:1).*parameters_of_rotor.pole_pairs_of_rotor;
            %% calculate the flux linkage in SD model
            Np = parameters_of_stator.number_of_phase;
            t = parameters_of_other_part.period_of_machine;
            if (mod(Ns, Np*t)==0 && Ps<Ns && Ps>0)
                parameters_of_performance.flux_linkage = Flux_linkage_calculation_based_on_linear_subdomain(parameters_of_rotor, parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics);
            else
                parameters_of_performance.flux_linkage =0;
            end
            %% record the data 
            n = n+1;
            SD.motor(n) =  Record_parameters_of_motor_subdomain(parameters_of_rotor, parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, parameters_of_performance);
        end
    end
end
toc

%% Data Processing
% transfer the struct of parameters into a table
table = struct2table(SD.motor);
% use the core parameters of the table;
SD_table = table(:,{'number_of_slot', 'pole_pairs_of_stator', 'pole_pairs_of_rotor', 'ratio_of_tooth', 'flux_linkage'});
% save the table in the mat docoment
writetable(SD_table, 'SD_table.txt');

