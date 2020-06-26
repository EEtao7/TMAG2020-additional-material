function current_of_every_phase= Current_of_every_phase( parameters_of_stator, parameters_of_rotor, parameters_of_time_harmonics, thetam_shift)
% This function is used to generate the current of every phase based on
% steady state analysis
%% Define some useful Parameters
Pr = parameters_of_rotor.pole_pairs_of_rotor;
Ns = parameters_of_stator.number_of_slot;
Np = parameters_of_stator.number_of_phase;
q = gcd(Pr,Ns);
current_amplitude = parameters_of_stator.current_amplitude;
time_harmonics_vector = parameters_of_time_harmonics.current_time_harmonics;

%% Calculate the current in every phase
current_of_every_phase = zeros(1, Np);
for m = 1: 1: Np
    n = 0;
    for time_harmonics = time_harmonics_vector
        n = n+1;
        if (mod(round(Pr/q), Np) == 1)
            % for the positive phase sequence
            current_of_every_phase(m) = current_of_every_phase(m)+current_amplitude(n)*cos(thetam_shift*time_harmonics-(m-1)/Np*2*pi);
        elseif(mod(round(Pr/q), Np) == (Np-1))
            % for the negative phase sequence
            current_of_every_phase(m) = current_of_every_phase(m)+current_amplitude(n)*cos(thetam_shift*time_harmonics+(m-1)/Np*2*pi);
        end
    end
end
end

