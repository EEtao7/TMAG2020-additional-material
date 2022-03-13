function Az = Az_subdomain_time_harmonics(radius, angle, parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, integrated_winding_matrix, MatrixCNmn, thetam, unit_coefficient)
% MatrixR is used to generate the magnetic vector potential Az in the motor
%% Define some useful parameters
Ns = parameters_of_stator.number_of_slot;
Pr = parameters_of_rotor.pole_pairs_of_rotor;
unit_coefficient = parameters_of_other_part.unit_coefficient;
r1 = parameters_of_rotor.radius_of_rotor/unit_coefficient;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/unit_coefficient;
r3 = parameters_of_stator.radius_of_stator/unit_coefficient;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/unit_coefficient;
q = parameters_of_other_part.period_of_machine;
delta3 = 2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics_vector = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
slot_max = parameters_of_stator.number_of_slot/q;
unit_slot = parameters_of_stator.unit_slots;
[number_of_phase, ~] = size(integrated_winding_matrix);

time_harmonics_vertor = parameters_of_time_harmonics.time_harmonics_collection;
time_harmonics_vertor = union(time_harmonics_vertor, -time_harmonics_vertor);

r = radius/unit_coefficient;
angle = mod(angle, 2*pi/q);

%% Compute the magnetic vector potential Az in the motor
% Initial some parameters
Az = 0;
angle_range_flag = 0;
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

m = 0;
for time_harmonics = time_harmonics_vertor
    m = m+1;
    %% Calculate the current density
    Jex = zeros(unit_slot, lamda_max+1);
    beta = 0;
    % When calculating the Az, we need to get the positive and negative time harmonics in the same time;
    current_time_harmonics_vector = parameters_of_time_harmonics.current_time_harmonics;
    if(ismember(abs(time_harmonics), current_time_harmonics_vector) == 1)
        [~, col_Jex] = find(current_time_harmonics_vector == abs(time_harmonics));
        Jex_amplitude = parameters_of_stator.Jex_amplitude(col_Jex);
        for slot = 1: 1: unit_slot
            % Decide the phase of current in the windings
            [row, col] = find(integrated_winding_matrix== slot | integrated_winding_matrix == -slot); % Decide the phase of the coils in every slots
            for number_of_layers = 1:1:2 % 2代表双层绕组的绕组层数
                for lamda = lamda_min: dlamda: lamda_max
                    n = round((lamda-lamda_min)/dlamda+1);
                    if(lamda == 0)% Decide the space harmonics lamda
                        if (mod(round(Pr/q), number_of_phase) == 1)
                            % for the positive phase sequence
                            alpha = (row(number_of_layers)-1)/number_of_phase*2*pi;
                        elseif(mod(round(Pr/q), number_of_phase) == (number_of_phase-1))
                            % for the negative phase sequence
                            alpha = -(row(number_of_layers)-1)/number_of_phase*2*pi;
                        end
                        Jex(slot, n) = Jex(slot, n)+0.5*Jex_amplitude*exp(sign(time_harmonics)*1i*(alpha-beta))*sign(integrated_winding_matrix(row(number_of_layers), col(number_of_layers))); % Decide the direction of current
                    else
                        Jex(slot, n) = 0;
                    end
                end
            end
        end
    end
    
    %% Calculate the Az in different region
    % Divide the region in the motor
    if(r>=r1 && r<r2)
        % for the PM region
        n = 0;
        for space_harmonics = space_harmonics_vector
            n = n+1;
            m_offset1 = 0;
            n1 = m_offset1+2*n-1;
            n2 = m_offset1+2*n;
            A1mn = MatrixCNmn(m, n1)*r^(-space_harmonics*q)+MatrixCNmn(m, n2)*r^(space_harmonics*q)+...
                P1n(Br_radial_megnetization_time_harmonics(parameters_of_rotor, parameters_of_time_harmonics, space_harmonics*q, time_harmonics), space_harmonics*q, r, thetam);
            Az = Az+A1mn*exp(1i*space_harmonics*q*angle)*exp(-1i*time_harmonics*thetam)*unit_coefficient;
        end
        
    elseif(r>=r2 && r<r3)
        % for the air-gap region
        n = 0;
        for space_harmonics = space_harmonics_vector
            n = n+1;
            m_offset2 = 2*(size_of_space_harmonics_vector);
            n1 = m_offset2+2*n-1;
            n2 = m_offset2+2*n;
            A2mn = MatrixCNmn(m, n1)*r^(-space_harmonics*q)+MatrixCNmn(m, n2)*r^(space_harmonics*q);
            Az = Az+A2mn*exp(1i*space_harmonics*q*angle)*exp(-1i*time_harmonics*thetam)*unit_coefficient;
        end
        
    elseif(r>=r3 && r<=r4)
        % for the slots region
        for slot = slot_min: dslot: slot_max
            alphaj = (slot-1)*2*pi/Ns+(pi/Ns);
            if(angle>alphaj-0.5*delta3 && angle<alphaj+0.5*delta3)
                for lamda = lamda_min: dlamda: lamda_max
                    k = 0;
                    while(slot>k*unit_slot)
                        k=k+1;
                    end
                    core_slot = mod(slot, unit_slot);
                    if(core_slot == 0)
                        core_slot = core_slot + unit_slot;
                    end
                    n = round((lamda-lamda_min)/dlamda+1);
                    m_offset3 = 4*(size_of_space_harmonics_vector);
                    n1 = m_offset3+(core_slot-1)*2*(lamda_max+1)+2*n-1;
                    n2 = m_offset3+(core_slot-1)*2*(lamda_max+1)+2*n;
                    if(lamda == 0)
                        A3mlamda = (MatrixCNmn(m, n1)+MatrixCNmn(m, n2)*log(r)+P2lamdaj(Jex(core_slot, n), lamda, r))*exp(1i*time_harmonics*2*pi/Ns*unit_slot*(k-1));
                    elseif(lamda~=0)
                        A3mlamda = MatrixCNmn(m, n1)*r^(-lamda*pi/delta3)+MatrixCNmn(m, n2)*r^(lamda*pi/delta3)+P2lamdaj(Jex(core_slot, n), lamda, r)*exp(1i*time_harmonics*2*pi/Ns*unit_slot*(k-1));
                    end
                    Az = Az+A3mlamda*cos((lamda*pi/delta3)*(angle-alphaj+0.5*delta3))*exp(-1i*time_harmonics*thetam)*unit_coefficient;
                end
                angle_range_flag = 1;
            end
        end
        if(angle_range_flag ~=1)
            disp('errors: the point is in the iron of the stator');
        end
    else
        % for the other region
        disp('errors: the point beyond the solve area');
    end
end

%% define some useful functions
    function p1n = P1n(Brn, harmonics_q, r, thetam)
        if(harmonics_q==1 || harmonics_q==-1)
            p1n = 1i*harmonics_q*Brn*r*log(r)*0.5*exp(-1i*harmonics_q*thetam);
        else
            p1n = 1i*harmonics_q*Brn*r/(1-harmonics_q^2)*exp(-1i*harmonics_q*thetam);
        end
    end

    function p2lamdaj = P2lamdaj(Jexlamdaj, lamda, r)
        mu0 = 4*pi*10^-7;
        constant = lamda*pi/delta3;
        if(constant==2 || constant==-2)
            p2lamdaj = -0.25*Jexlamdaj*r^2*mu0*log(r);
        else
            p2lamdaj = Jexlamdaj*r^2*mu0/(constant^2-4);
        end
    end

end