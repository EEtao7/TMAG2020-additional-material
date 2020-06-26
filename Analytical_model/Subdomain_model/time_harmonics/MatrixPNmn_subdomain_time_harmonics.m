function matrixPNmn = MatrixPNmn_subdomain_time_harmonics(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics, time_harmonics)
% MatrixPNmn_subdomain_time_harmonics is used to generate the matrixPNmn
%% Define some useful parameters
mur = parameters_of_rotor.mur;
Ns = parameters_of_stator.number_of_slot;
unit_coefficient = parameters_of_other_part.unit_coefficient;
r1 = parameters_of_rotor.radius_of_rotor/unit_coefficient;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/unit_coefficient;
r3 = parameters_of_stator.radius_of_stator/unit_coefficient;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/unit_coefficient;
r_ref1 = parameters_of_other_part.r_reference1/unit_coefficient;
r_ref2 = parameters_of_other_part.r_reference2/unit_coefficient;
q = parameters_of_other_part.period_of_machine;
delta3 = 2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
unit_slot = parameters_of_stator.unit_slots;
N_total = 4*(size_of_space_harmonics)+2*unit_slot*(lamda_max+1);

%% Generate the matrixPN
% Initial some parameters
matrixPNmn = zeros(N_total, N_total);
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

%1. for boundary condition r = r1;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset1 = 0;
    m = m_offset1+n;
    matrixPNmn(m, :) = [zeros(1,(n-1)*2), -(space_harmonics*q)*(r1/r_ref1)^(-(space_harmonics*q)-1)/r_ref1, (space_harmonics*q)*(r1/r_ref1)^((space_harmonics*q)-1)/r_ref1, zeros(1, N_total-2*n)];
end

%2. for boundary condition1 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset2 = size_of_space_harmonics;
    m = m_offset2+n;
    matrixPNmn(m, :) = [zeros(1,(n-1)*2), -(space_harmonics*q)*(r2/r_ref1)^(-(space_harmonics*q)-1)/r_ref1, (space_harmonics*q)*(r2/r_ref1)^((space_harmonics*q)-1)/r_ref1, zeros(1, 2*(size_of_space_harmonics)-2*n),...
        zeros(1,(n-1)*2), mur*(space_harmonics*q)*(r2/r_ref1)^(-(space_harmonics*q)-1)/r_ref1, -mur*(space_harmonics*q)*(r2/r_ref1)^(space_harmonics*q-1)/r_ref1, zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
end

%3. for boundary condition2 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset3 = 2*(size_of_space_harmonics);
    m = m_offset3+n;
    matrixPNmn(m, :) = [zeros(1,(n-1)*2), (r2/r_ref1)^(-(space_harmonics*q)), (r2/r_ref1)^(space_harmonics*q), zeros(1, 2*(size_of_space_harmonics)-2*n),...
        zeros(1,(n-1)*2), -(r2/r_ref1)^(-(space_harmonics*q)), -(r2/r_ref1)^(space_harmonics*q), zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
end

%4. for the boudary condition r = r4;
for slot = slot_min: dslot: unit_slot
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset4 = 3*(size_of_space_harmonics);
        m = m_offset4+(slot-1)*(lamda_max+1)+n;
        if(lamda ~= 0)
            matrixPNmn(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -(lamda*pi/delta3)*(r4/r_ref2)^(-(lamda*pi/delta3)-1)/r_ref2, (lamda*pi/delta3)*(r4/r_ref2)^((lamda*pi/delta3)-1)/r_ref2,...
                zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
        elseif(lamda == 0)
            matrixPNmn(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                0, 1/r4, zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
        end
    end
end

%5. for the boudary condition1 r = r3;
for slot = slot_min: dslot: unit_slot
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset5 = 3*(size_of_space_harmonics)+unit_slot*(lamda_max+1);
        m = m_offset5+(slot-1)*(lamda_max+1)+n;
        if(lamda ~= 0)
            matrixPNmn(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -(r3/r_ref2)^(-(lamda*pi/delta3)), -(r3/r_ref2)^(lamda*pi/delta3), zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
            h = 0;
            for space_harmonics = space_harmonics_vector
                h = h+1;
                matrixPNmn(m, :) = matrixPNmn(m, :)+[zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(h-1)),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*(r3/r_ref1)^(-space_harmonics*q),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*(r3/r_ref1)^(space_harmonics*q),...
                    zeros(1, N_total-2*(size_of_space_harmonics)-2*h)];
            end
        elseif(lamda == 0)
            matrixPNmn(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -1, -log(r3), zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
            h = 0;
            for space_harmonics = space_harmonics_vector
                h = h+1;
                matrixPNmn(m, :) = matrixPNmn(m, :)+[zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(h-1)),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*(r3/r_ref1)^(-space_harmonics*q),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*(r3/r_ref1)^(space_harmonics*q),...
                    zeros(1, N_total-2*(size_of_space_harmonics)-2*h)];
            end
        end
    end
end

%5. for the boudary condition2 r = r3;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset6 = 3*(size_of_space_harmonics)+2*unit_slot*(lamda_max+1);
    m = m_offset6+n;
    matrixPNmn(m, :) = [zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(n-1)),...
        (-space_harmonics*q)*(r3/r_ref1)^(-space_harmonics*q-1)/r_ref1, (space_harmonics*q)*(r3/r_ref1)^(space_harmonics*q-1)/r_ref1,...
        zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
    offset_coefficient = 0;% initial the offset coefficient
    for k = 0: 1: (Ns/(unit_slot*q)-1)
        offset_coefficient =  offset_coefficient+exp(1i*(time_harmonics-space_harmonics*q)*2*pi/(Ns/unit_slot)*k);
    end
    if(abs(offset_coefficient)>1e-10)% if offset coefficient is close to 0, we don't use this small value
        for slot = slot_min: dslot: unit_slot
            for lamda = lamda_min: dlamda: lamda_max
                l = round((lamda-lamda_min)/dlamda+1);
                if(lamda ~= 0)
                    matrixPNmn(m, :) = matrixPNmn(m, :)+[zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(l-1)),...
                        offset_coefficient*(-q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*(-lamda*pi/delta3)*(r3/r_ref2)^(-lamda*pi/delta3-1))/r_ref2,...
                        offset_coefficient*(-q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*(lamda*pi/delta3)*(r3/r_ref2)^(lamda*pi/delta3-1))/r_ref2,...
                        zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*l)];
                elseif(lamda == 0)
                    matrixPNmn(m, :) = matrixPNmn(m, :)+[zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(l-1)),...
                        0, offset_coefficient*(-q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*1/r3),...
                        zeros(1, N_total-4*(size_of_space_harmonics)- 2*(slot-1)*(lamda_max+1)-2*l)];
                end
            end
        end
    end
end

%% Assistant function
    function f2 = f2(lamda)
        if(lamda==0)
            f2 = 0.5;
        elseif(lamda~=0)
            f2 = 1;
        end
    end
    function m1 = M1(harmonics, lamda, j)
        alphaj = (j-1)*2*pi/Ns+pi/Ns;
        m1 = F(delta3*q*harmonics-pi*lamda)*exp(1i*(0.5*pi*lamda-alphaj*q*harmonics))*delta3;
    end
    function m2 = M2(harmonics, lamda, j)
        alphaj = (j-1)*2*pi/Ns+pi/Ns;
        m2 = F(delta3*q*harmonics-pi*lamda)*exp(1i*(alphaj*q*harmonics-0.5*pi*lamda));
    end
    function F = F(x)
        if(x==0)
            F = 1;
        elseif(x~=0)
            F = sin(0.5*x)/(0.5*x);
        end
    end

end

