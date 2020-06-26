%% Flux linkage Test Based on Linear FEM  (Vernier PMSM)
% This program is designed to test the flux linkage in the Vernier PMSM
% based on the linear FEM
% Process Oriented Programming
% EEtao
% tangchentao@zju.edu.cn
%
clc; clear;

%% Induce the functions and packages
addpath('.\Motor_Objects')
addpath('.\Motor_Component');
addpath('.\Design_Engineer');
addpath('.\Test_Engineer');
addpath('.\Material');
addpath('.\Motor_design_toolbox');
addpath('.\FFT_analysis');
import Winding_Matrix_Package.*;
import Air_gap_Package.*;
import Rotor_Package.*;
import Stator_Package.*;
import Windings_Package.*;

%% Initial Parameters of motor (SI)
Nonlinear_flag = false; % true
    % when flag = true, FEM use the nonlinear material 
    % when flag = false, FEM use the linear material
    
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
parameters_of_stator.depth_of_slot = 10e-3;
parameters_of_stator.pole_pairs_of_stator = 5;
parameters_of_stator.number_of_phase = 3;
parameters_of_stator.turns_of_phase = 200;
parameters_of_stator.Jex_amplitude = 0; %(A/m^2)

% parameters of other part
parameters_of_other_part.length_of_stack = 10e-3;

%% Some useful coefficient
unit_coefficient = 1e3;

%% Transfer the parameters to FEM model
%% Design the detail parameters of  Stator
Outer_Auxiliary_Tooth_Stator = OuterAuxiliaryToothStator();
% basic parameters
Outer_Auxiliary_Tooth_Stator.name = 'Outer_Auxiliary_Stator_';
Outer_Auxiliary_Tooth_Stator.group_number_of_stator = 2;
Outer_Auxiliary_Tooth_Stator.precision_of_stator_iron = 5;
Outer_Auxiliary_Tooth_Stator.outer_radius = 100;
% special parameters of SufacemountedPM_Rotor
Outer_Auxiliary_Tooth_Stator.material_of_Iron = Material_Iron();
Outer_Auxiliary_Tooth_Stator.nonlinear_flag = Nonlinear_flag;

% adjust parameters
Outer_Auxiliary_Tooth_Stator.pole_pairs_of_stator = parameters_of_stator.pole_pairs_of_stator ;
Outer_Auxiliary_Tooth_Stator.number_of_slots  = parameters_of_stator.number_of_slot;
Outer_Auxiliary_Tooth_Stator.number_of_auxiliary_tooth  = 1;
Outer_Auxiliary_Tooth_Stator.ratio_of_tooth = 0.2;
Outer_Auxiliary_Tooth_Stator.ratio_of_auxiliary_tooth = parameters_of_stator.ratio_of_tooth;
Outer_Auxiliary_Tooth_Stator.slot_depth = 19.5;

%% Design the detail parameters of Windings
Outer_Auxiliary_Tooth_Stator.auxiliary_slot_depth = 1.5;
Outer_Auxiliary_Tooth_Stator.thick_of_auxiliary_yoke_stator = 2;
Outer_Auxiliary_Tooth_Stator.adjust_coefficient_of_thickness_of_auxiliary_yoke_stator = 0.5;
Outer_Auxiliary_Tooth_Stator.thick_of_yoke_stator = 20;

%% parameters_of_windingsdings
Windings = GeneralWindings();
% basic parameters
Windings.name =  'GeneralWindings_';
Windings.precision_of_windings = 4;
Windings.group_number_of_winding = 10;
Windings.group_number_of_windingA = 11;
Windings.group_number_of_windingB = 12;
Windings.group_number_of_windingC = 13;
Windings.current_of_windings = Current_of_Windings();
Windings.stator = Outer_Auxiliary_Tooth_Stator;
Windings.pole_pairs_of_stator = Windings.stator.pole_pairs_of_stator;
Windings.number_of_slots =Windings.stator.number_of_slots;
Windings.number_of_phase = parameters_of_stator.number_of_phase;
% special parameters of Fractional Slots Windings
Windings.series_flag = 1;
Windings.number_of_layers =2;
Windings.turns_of_phase = parameters_of_stator.turns_of_phase;
unit_windingmatrix_of_phase_A= Unit_WindingMatrix_of_PhaseA(Windings);
Windings.Winding_Matrix_PhaseA = unit_windingmatrix_of_phase_A;


%% Design the detail parameters of AirGap
Inner_AirGap = AGE_Inner_Rotor();
% basic parameters
Inner_AirGap.name =  'AGE_Inner_Rotor_';
Inner_AirGap.group_number_of_AirGap = 3;
Inner_AirGap.precision_of_AirGap = 1;
% special parameters of AGE
Inner_AirGap. stator= Outer_Auxiliary_Tooth_Stator;
Inner_AirGap.outer_AGE_angle = 0;
Inner_AirGap.inner_AGE_angle = 0;
% adjust parameters
Inner_AirGap.air_gap_length = 1;
Inner_AirGap.AGE_flag = 6;

%% Design the detail parameters of Rotor
% SurfacemountedPM rotor
Inner_SP_Rotor = SurfacemountedPM_Inner_Rotor();
% basic parameters
Inner_SP_Rotor.name =  'SurfacemountedPM_Inner_Rotor_';
Inner_SP_Rotor.group_number_of_rotor =0;
Inner_SP_Rotor.precision_of_rotor_iron = 5;
% special parameters of SufacemountedPM_Rotor
Inner_SP_Rotor.precision_of_PM = 2;
Inner_SP_Rotor.group_number_of_PM = 1;
Inner_SP_Rotor.material_of_PM = Material_PM();
Inner_SP_Rotor.material_of_Iron = Material_Iron();
Inner_SP_Rotor.nonlinear_flag = Nonlinear_flag;

% adjust parameters
Inner_SP_Rotor.pole_of_pairs_Rotor = parameters_of_rotor.pole_pairs_of_rotor ;
Inner_SP_Rotor.radius_of_PM = Outer_Auxiliary_Tooth_Stator.inner_radius-Inner_AirGap.air_gap_length;
Inner_SP_Rotor.thickness_of_PM = parameters_of_rotor.thickness_of_PM*unit_coefficient;
Inner_SP_Rotor.inner_radious = 20;
Inner_SP_Rotor.ratio_of_rotor_PM_bottom = parameters_of_rotor.ratio_of_magnet;
Inner_SP_Rotor.ratio_of_rotor_PM_top = 1.0;

%% Design a PMSM in FEM
PMSM_Flux_linkage_Test = Motor();
%define every part of PMSM
PMSM_Flux_linkage_Test.stator = Outer_Auxiliary_Tooth_Stator;
PMSM_Flux_linkage_Test.windings = Windings;
PMSM_Flux_linkage_Test.air_gap = Inner_AirGap;
PMSM_Flux_linkage_Test.rotor = Inner_SP_Rotor;
PMSM_Flux_linkage_Test.stack_length = parameters_of_other_part.length_of_stack*unit_coefficient;

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

% Calculate and save the data in the FEM model
tic
n = 0;
for Ns = Ns_range
    parameters_of_stator.number_of_slot = Ns;
    Outer_Auxiliary_Tooth_Stator.number_of_slots  = parameters_of_stator.number_of_slot; % change the parameters in the FEM
    for Pr = Pr_range
        parameters_of_rotor.pole_pairs_of_rotor = Pr;
        Inner_SP_Rotor.pole_of_pairs_Rotor = parameters_of_rotor.pole_pairs_of_rotor ; % change the parameters in the FEM
        Ps = Ns-Pr; % for the Vernier PMSM
        parameters_of_stator.pole_pairs_of_stator = Ps;
        Outer_Auxiliary_Tooth_Stator.pole_pairs_of_stator = parameters_of_stator.pole_pairs_of_stator ; % change the parameters in the FEM
        for Ratio_of_tooth = Ratio_of_tooth_range
            parameters_of_stator.ratio_of_tooth = Ratio_of_tooth;
            Outer_Auxiliary_Tooth_Stator.ratio_of_auxiliary_tooth = parameters_of_stator.ratio_of_tooth; % change the parameters in the FEM
             %% adjust some assistant parameters with the key parameters
            Windings.pole_pairs_of_stator = Windings.stator.pole_pairs_of_stator;
            Windings.number_of_slots = Windings.stator.number_of_slots;
            unit_windingmatrix_of_phase_A= Unit_WindingMatrix_of_PhaseA(Windings);
            Windings.Winding_Matrix_PhaseA = unit_windingmatrix_of_phase_A;
            %% calculate the flux linkage in FEM model
            Np = parameters_of_stator.number_of_phase;
            parameters_of_other_part.period_of_machine = gcd(parameters_of_rotor.pole_pairs_of_rotor, parameters_of_stator.number_of_slot);
            t = parameters_of_other_part.period_of_machine;
            if (mod(Ns, Np*t)==0 && Ps<Ns)
                parameters_of_performance.flux_linkage = Flux_linkage_calculation_based_on_linear_FEM(PMSM_Flux_linkage_Test);
            else
                parameters_of_performance.flux_linkage =0;
            end
            %% record the data 
            n = n+1;
            FEM.motor(n) =  Record_parameters_of_motor_FEM(parameters_of_rotor, parameters_of_stator, parameters_of_other_part,  parameters_of_performance);
        end
    end
end
toc

%% Data Processing
% transfer the struct of parameters into a table
table = struct2table(FEM.motor);
% use the core parameters of the table;
FEM_table = table(:,{'number_of_slot', 'pole_pairs_of_stator', 'pole_pairs_of_rotor', 'ratio_of_tooth', 'flux_linkage'});
% save the table in the txt docoment
writetable(FEM_table, 'FEM_table.txt');


