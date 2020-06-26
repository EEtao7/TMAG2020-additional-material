function matrixRN = MatrixRN_subdomain_time_stepping(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, integrated_winding_matrix, matrixJ_of_every_phase, thetam)
% This function is used to generate the matrixRN
%% Define some useful parameters
Ns = parameters_of_stator.number_of_slot;
r1 = parameters_of_rotor.radius_of_rotor/parameters_of_other_part.r_reference;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/parameters_of_other_part.r_reference;
r3 = parameters_of_stator.radius_of_stator/parameters_of_other_part.r_reference;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/parameters_of_other_part.r_reference;
q = parameters_of_other_part.period_of_machine;
delta3 = parameters_of_stator.radius_of_stator*2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot/parameters_of_other_part.r_reference;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics_vector = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
slot_max = parameters_of_stator.number_of_slot/q;
N_total = 4*(size_of_space_harmonics_vector)+2*slot_max*(lamda_max+1);

[number_of_phase, ~] = size(integrated_winding_matrix); % for the jex

%% Generate the matrixRN
% Initial some parameters
matrixRN = zeros(N_total, 1);
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
                if(lamda == 0)% Decide the space harmonics lamda
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

%% Boundary conditions in the matrixRN

% 1. for boundary condition r = r1;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset1 = 0;
    m = m_offset1+n;
   matrixRN(m,1) = -D1n(Br_radial_megnetization(parameters_of_rotor, space_harmonics*q), space_harmonics*q, r1, thetam);
end

% 2. for boundary condition1 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset2 = size_of_space_harmonics_vector;
    m = m_offset2+n;
    matrixRN(m,1) = -D1n(Br_radial_megnetization(parameters_of_rotor, space_harmonics*q), space_harmonics*q, r2, thetam);
end

% 3. for boundary condition2 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset3 = 2*(size_of_space_harmonics_vector);
    m = m_offset3+n;
    matrixRN(m,1) = -P1n(Br_radial_megnetization(parameters_of_rotor, space_harmonics*q), space_harmonics*q, r2, thetam);
end

% 4. for the boudary condition r = r4;
for slot = slot_min: dslot: slot_max
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset4 = 3*(size_of_space_harmonics_vector);
        m = m_offset4+(slot-1)*(lamda_max+1)+n;    
        matrixRN(m,1) = -D2lamdaj(Jex(slot, n), lamda, r4);
    end
end

% 5. for the boudary condition1 r = r3;
for slot = slot_min: dslot: slot_max
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset5 = 3*(size_of_space_harmonics_vector)+slot_max*(lamda_max+1);
        m = m_offset5+(slot-1)*(lamda_max+1)+n;
        matrixRN(m,1) = P2lamdaj(Jex(slot,n), lamda, r3);
    end
end

% 6. for the boudary condition2 r = r3;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset6 = 3*(size_of_space_harmonics_vector)+2*slot_max*(lamda_max+1);
    m = m_offset6+n;
    for slot = slot_min: dslot: slot_max
        for lamda = lamda_min: dlamda: lamda_max
            l = round((lamda-lamda_min)/dlamda+1);
            matrixRN(m,1) = matrixRN(m,1)+q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*D2lamdaj(Jex(slot, l), lamda, r3);
        end
    end
end

%% Define some useful functions
    function p1n = P1n(Brn, harmonics_q, r, thetam)
        if(harmonics_q==1 || harmonics_q==-1)
            p1n = 1i*harmonics_q*Brn*r*log(r)*0.5*exp(-1i*harmonics_q*thetam);
        else
            p1n = 1i*harmonics_q*Brn*r/(1-harmonics_q^2)*exp(-1i*harmonics_q*thetam);
        end
    end

    function d1n = D1n(Brn, harmonics_q, r, thetam)
        if(harmonics_q==1 || harmonics_q==-1)
            d1n = 1i*harmonics_q*Brn*(log(r)+1)*0.5*exp(-1i*harmonics_q*thetam);
        else
            d1n = 1i*harmonics_q*Brn*1/(1-harmonics_q^2)*exp(-1i*harmonics_q*thetam);
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

    function d2lamdaj = D2lamdaj(Jexlamdaj, lamda, r)
        mu0 = 4*pi*10^-7;
        k = lamda*pi/delta3;
        if(k==2 || k==-2)
            d2lamdaj = -0.25*Jexlamdaj*mu0*(2*r*log(r)+r);
        else
            d2lamdaj = Jexlamdaj*2*r*mu0/(k^2-4);
        end
    end

    function m1 = M1(harmonics, lamda, j)
        alphaj = (j-1)*2*pi/Ns+pi/Ns;
        m1 = F(delta3*q*harmonics-pi*lamda)*exp(1i*(0.5*pi*lamda-alphaj*q*harmonics))*delta3;
    end

    function F = F(x)
        if(x==0)
            F = 1;
        elseif(x~=0)
            F = sin(0.5*x)/(0.5*x);
        end
    end

end

