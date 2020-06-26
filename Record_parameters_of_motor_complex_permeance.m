function motor = Record_parameters_of_motor_complex_permeance(parameters_of_rotor, parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, parameters_of_performance)
% This function is used to initialize the parameters of motor in the complex
% permeance model

%% Record the parameters of motor
% parameters of rotor
motor.B_remanence = parameters_of_rotor.B_remanence;
motor.mur = parameters_of_rotor.mur;
motor.pole_pairs_of_rotor = parameters_of_rotor.pole_pairs_of_rotor;
motor.radius_of_rotor = parameters_of_rotor.radius_of_rotor;
motor.ratio_of_magnet = parameters_of_rotor.ratio_of_magnet;
motor.thickness_of_PM = parameters_of_rotor.thickness_of_PM;

% parameters of stator
motor.radius_of_stator  = parameters_of_stator.radius_of_stator;
motor.ratio_of_tooth  = parameters_of_stator.ratio_of_tooth;
motor.number_of_slot = parameters_of_stator.number_of_slot;
motor.pole_pairs_of_stator = parameters_of_stator.pole_pairs_of_stator;
motor.number_of_phase = parameters_of_stator.number_of_phase;
motor.turns_of_phase = parameters_of_stator.turns_of_phase;
motor.current_amplitude = parameters_of_stator.current_amplitude; 

% parameters of other part
motor.length_of_stack = parameters_of_other_part.length_of_stack;
motor.period_of_machine = parameters_of_other_part.period_of_machine;
motor.air_gap = parameters_of_other_part.air_gap;
motor.radius = parameters_of_other_part.radius; 

% parameters of time harmonics
motor.current_time_harmonics = parameters_of_time_harmonics.current_time_harmonics ;
motor.PM_time_harmonics = parameters_of_time_harmonics.PM_time_harmonics ;

% parameters of space harmonics
motor.space_harmonics_constant_max = parameters_of_space_harmonics.space_harmonics_constant_max;

% parameters of performance
motor.flux_linkage = parameters_of_performance.flux_linkage;

end

