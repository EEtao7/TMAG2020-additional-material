function Az = Az_subdomain_time_stepping(radius, angle, parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, matrixCN, integrated_winding_matrix, matrixJ_of_every_phase, thetam)
% MatrixR is used to generate the magnetic vector potential Az in the motor
%% Define some useful parameters
Ns = parameters_of_stator.number_of_slot;
r1 = parameters_of_rotor.radius_of_rotor/parameters_of_other_part.r_reference;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/parameters_of_other_part.r_reference;
r3 = parameters_of_stator.radius_of_stator/parameters_of_other_part.r_reference;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/parameters_of_other_part.r_reference;
q = parameters_of_other_part.period_of_machine;
delta3 = 2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics_vector = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
slot_max = parameters_of_stator.number_of_slot/q;

r = radius/parameters_of_other_part.r_reference;
angle = mod(angle, 2*pi/q);
r_reference = parameters_of_other_part.r_reference;

[number_of_phase, ~] = size(integrated_winding_matrix); % for the jex

%% Compute the magnetic vector potential Az in the motor
% Initial some parameters
Az = 0;
angle_range_flag = 0;
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

%% Generate the current density matrix for the windings
Jex = zeros( slot_max, lamda_max+1);
% check the matrix J of every phase
if(number_of_phase == length(matrixJ_of_every_phase))
    for slot = 1: 1: slot_max
        [row, col] = find(integrated_winding_matrix== slot | integrated_winding_matrix == -slot); % Decide the phase of the coils in every slots
        for number_of_layer = 1:1:2 % 2代表双层绕组的绕组层数
            for lamda = lamda_min: dlamda: lamda_max
                n = round((lamda-lamda_min)/dlamda+1);
                if(lamda == 0) % Decide the space harmonics lamda
                    Jex(slot, n) = Jex(slot, n)+matrixJ_of_every_phase(row(number_of_layer))*sign(integrated_winding_matrix(row(number_of_layer), col(number_of_layer))); % Decide the direction of current
                else
                    Jex(slot, n) = 0;
                end
            end
        end
    end
else
    disp('Errors in the J matrix of every phase');
end

%% Calculate the Az in the different regions
% Divide the region in the motor
if(r>=r1 && r<r2)
    % for the PM region
    n = 0;
    for space_harmonics = space_harmonics_vector
        n = n+1;
        m_offset1 = 0;
        m1 = m_offset1+2*n-1;
        m2 = m_offset1+2*n;
        A1n = matrixCN(m1, 1)*r^(-space_harmonics*q)+matrixCN(m2, 1)*r^(space_harmonics*q)+P1n(Br_radial_megnetization(parameters_of_rotor, space_harmonics*q), space_harmonics*q, r, thetam);
        Az = Az+A1n*exp(1i*space_harmonics*q*angle)*r_reference;
    end
    
elseif(r>=r2 && r<r3)
    % for the Air-gap region
    n = 0;
    for space_harmonics = space_harmonics_vector
        n = n+1;
        m_offset2 = 2*(size_of_space_harmonics_vector);
        m1 = m_offset2+2*n-1;
        m2 = m_offset2+2*n;
        A2n = matrixCN(m1, 1)*r^(-space_harmonics*q)+matrixCN(m2, 1)*r^(space_harmonics*q);
        Az = Az+A2n*exp(1i*space_harmonics*q*angle)*r_reference;
    end
    
elseif(r>=r3 && r<=r4)
    % for the Slots region
    for j = slot_min: dslot: slot_max
        alphaj = (j-1)*2*pi/Ns+(pi/Ns);
        if(angle>alphaj-0.5*delta3 && angle<alphaj+0.5*delta3)
            for lamda = lamda_min: dlamda: lamda_max
                n = round((lamda-lamda_min)/dlamda+1);
                m_offset3 = 4*(size_of_space_harmonics_vector);
                m1 = m_offset3+(j-1)*2*(lamda_max+1)+2*n-1;
                m2 = m_offset3+(j-1)*2*(lamda_max+1)+2*n;
                if(lamda == 0)
                    A3lamda = matrixCN(m1,1)+matrixCN(m2,1)*log(r)+P2lamdaj(Jex(j, n), lamda, r);
                elseif(lamda~=0)
                    A3lamda = matrixCN(m1,1)*r^(-lamda*pi/delta3)+matrixCN(m2,1)*r^(lamda*pi/delta3)+P2lamdaj(Jex(j, n), lamda, r);
                end
                Az = Az+A3lamda*cos((lamda*pi/delta3)*(angle-alphaj+0.5*delta3))*r_reference;
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

%% Define some useful functions
    function p1n = P1n(Brn, harmonics_q, r, thetam)
        if(harmonics_q==1 || harmonics_q==-1)
            p1n = 1i*harmonics_q*Brn*r*log(r)*0.5*exp(-1i*harmonics_q*thetam);
        else
            p1n = 1i*harmonics_q*Brn*r/(1-harmonics_q^2)*exp(-1i*harmonics_q*thetam);
        end
    end

    function p2lamdaj = P2lamdaj(Jexlamdaj, lamda, r)
        mu0 = 4*pi*10^-7;
        k = lamda*pi/delta3;
        if(k==2 || k==-2)
            p2lamdaj = -0.25*Jexlamdaj*r^2*mu0*log(r);
        else
            p2lamdaj = Jexlamdaj*r^2*mu0/(k^2-4);
        end
    end

end