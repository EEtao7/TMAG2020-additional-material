function matrixRNmn = MatrixRNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics, time_harmonics, integrated_winding_matrix)
% MatrixR is used to generate the matrixRN
%% Define some useful parameters
Ns = parameters_of_stator.number_of_slot;
Pr = parameters_of_rotor.pole_pairs_of_rotor;
unit_coefficient = parameters_of_other_part.unit_coefficient;
r1 = parameters_of_rotor.radius_of_rotor/unit_coefficient;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/unit_coefficient;
r3 = parameters_of_stator.radius_of_stator/unit_coefficient;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/unit_coefficient;
q = parameters_of_other_part.period_of_machine;
delta3 = parameters_of_stator.radius_of_stator*2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics_vector = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
unit_slot = parameters_of_stator.unit_slots;
N_total = 4*(size_of_space_harmonics_vector)+2*unit_slot*(lamda_max+1);

[number_of_phase, ~] = size(integrated_winding_matrix);

%% Generate the matrixRN
% Initial some parameters
matrixRNmn = zeros(N_total, 1);
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

%% Generate the current density matrix for the windings
Jex = zeros( unit_slot, lamda_max+1);
beta = 0;

% current_time_harmonics_vector = union(parameters_of_time_harmonics.current_time_harmonics, -1*parameters_of_time_harmonics.current_time_harmonics);
current_time_harmonics_vector = parameters_of_time_harmonics.current_time_harmonics; % just calculate the positive time harmonics
if(ismember(time_harmonics, current_time_harmonics_vector) == 1)
    [~, col_Jex] = find(current_time_harmonics_vector==time_harmonics);
    Jex_amplitude = parameters_of_stator.Jex_amplitude(col_Jex);
    for slot = 1: 1: unit_slot
        % Decide the phase of current in the windings
        [row, col] = find(integrated_winding_matrix== slot | integrated_winding_matrix == -slot); % Decide the phase of the coils in every slots
        for number_of_layer = 1:1:2 % 2代表双层绕组的绕组层数
            for lamda = lamda_min: dlamda: lamda_max
                n = round((lamda-lamda_min)/dlamda+1);
                if(lamda == 0)% Decide the space harmonics lamda
                    if (mod(round(Pr/q), number_of_phase) == 1)
                        % for the positive phase sequence
                        alpha = (row(number_of_layer)-1)/number_of_phase*2*pi;
                    elseif(mod(round(Pr/q), number_of_phase) == (number_of_phase-1))
                        % for the negative phase sequence
                        alpha = -(row(number_of_layer)-1)/number_of_phase*2*pi;
                    end
                    Jex(slot, n) = Jex(slot, n)+0.5*Jex_amplitude*exp(sign(time_harmonics)*1i*(alpha-beta))*sign(integrated_winding_matrix(row(number_of_layer), col(number_of_layer))); % Decide the direction of current
                else
                    Jex(slot, n) = 0;
                end
            end
        end
    end
end


%1. for boundary condition r = r1;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset1 = 0;
    m = m_offset1+n;
    matrixRNmn(m,1) = -D1mn(Br_radial_megnetization_time_harmonics(parameters_of_rotor, parameters_of_time_harmonics, space_harmonics*q, time_harmonics), space_harmonics*q, r1);
end

%2. for boundary condition1 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset2 = size_of_space_harmonics_vector;
    m = m_offset2+n;
    matrixRNmn(m,1) = -D1mn(Br_radial_megnetization_time_harmonics(parameters_of_rotor, parameters_of_time_harmonics, space_harmonics*q, time_harmonics), space_harmonics*q, r2);
end

%3. for boundary condition2 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset3 = 2*(size_of_space_harmonics_vector);
    m = m_offset3+n;
    matrixRNmn(m,1) = -P1mn(Br_radial_megnetization_time_harmonics(parameters_of_rotor, parameters_of_time_harmonics, space_harmonics*q, time_harmonics), space_harmonics*q, r2);
end

%4. for the boudary condition r = r4;
for slot = slot_min: dslot: unit_slot
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset4 = 3*(size_of_space_harmonics_vector);
        m = m_offset4+(slot-1)*(lamda_max+1)+n;
        matrixRNmn(m,1) = -D2mlamdaj(Jex(slot, n), lamda, r4);
    end
end

%5. for the boudary condition1 r = r3;
for slot = slot_min: dslot: unit_slot
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset5 = 3*(size_of_space_harmonics_vector)+unit_slot*(lamda_max+1);
        m = m_offset5+(slot-1)*(lamda_max+1)+n;
        matrixRNmn(m,1) = P2mlamdaj(Jex(slot,n), lamda, r3);
    end
end

n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset6 = 3*(size_of_space_harmonics_vector)+2*unit_slot*(lamda_max+1);
    m = m_offset6+n;
    offset_coefficient_current = 0;% initial the offset coefficient
    for s = 0: 1: (Ns/(unit_slot*q)-1)
        offset_coefficient_current =  offset_coefficient_current+exp(1i*(time_harmonics-space_harmonics*q)*2*pi/(Ns/unit_slot)*s);
    end
    if(abs(offset_coefficient_current)>1e-10)% if offset coefficient is close to 0, we don't use this small value
        for slot = slot_min: dslot: unit_slot
            for lamda = lamda_min: dlamda: lamda_max
                l = round((lamda-lamda_min)/dlamda+1);
                matrixRNmn(m,1) = matrixRNmn(m,1)+offset_coefficient_current*(q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*D2mlamdaj(Jex(slot, l), lamda, r3));
            end
        end
    end
end

%% define some useful functions
    function p1mn = P1mn(Brmn, space_harmonics_q, r)
        if(space_harmonics_q==1 || space_harmonics_q==-1)
            p1mn = 1i*space_harmonics_q*Brmn*r*log(r)*0.5;
        else
            p1mn = 1i*space_harmonics_q*Brmn*r/(1-space_harmonics_q^2);
        end
    end

    function d1mn = D1mn(Brmn, space_harmonics_q, r)
        if(space_harmonics_q==1 || space_harmonics_q==-1)
            d1mn = 1i*space_harmonics_q*Brmn*(log(r)+1)*0.5;
        else
            d1mn = 1i*space_harmonics_q*Brmn*1/(1-space_harmonics_q^2);
        end
    end

    function p2mlamdaj = P2mlamdaj(Jexlamdaj, lamda, r)
        mu0 = 4*pi*10^-7;
        k = lamda*pi/delta3;
        if(k==2 || k==-2)
            p2mlamdaj = -0.25*Jexlamdaj*r^2*mu0*log(r);
        else
            p2mlamdaj = Jexlamdaj*r^2*mu0/(k^2-4);
        end
    end

    function d2mlamdaj = D2mlamdaj(Jexlamdaj, lamda, r)
        mu0 = 4*pi*10^-7;
        k = lamda*pi/delta3;
        if(k==2 || k==-2)
            d2mlamdaj = -0.25*Jexlamdaj*mu0*(2*r*log(r)+r);
        else
            d2mlamdaj = Jexlamdaj*2*r*mu0/(k^2-4);
        end
    end

    function m1 = M1(space_harmonics, lamda, j)
        alphaj = (j-1)*2*pi/Ns+pi/Ns;
        m1 = F(delta3*q*space_harmonics-pi*lamda)*exp(1i*(0.5*pi*lamda-alphaj*q*space_harmonics))*delta3;
    end

    function F = F(x)
        if(x==0)
            F = 1;
        elseif(x~=0)
            F = sin(0.5*x)/(0.5*x);
        end
    end

end

