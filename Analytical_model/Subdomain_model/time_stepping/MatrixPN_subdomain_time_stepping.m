function matrixPN = MatrixPN_subdomain_time_stepping(parameters_of_stator, parameters_of_rotor, parameters_of_other_part, parameters_of_space_harmonics)
% This function is used to generate the matrixPN
%% Define some useful parameters
mur = parameters_of_rotor.mur;
Ns = parameters_of_stator.number_of_slot;
r1 = parameters_of_rotor.radius_of_rotor/parameters_of_other_part.r_reference;
r2 = (parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM)/parameters_of_other_part.r_reference;
r3 = parameters_of_stator.radius_of_stator/parameters_of_other_part.r_reference;
r4 = (parameters_of_stator.radius_of_stator+parameters_of_stator.depth_of_slot)/parameters_of_other_part.r_reference;
q = parameters_of_other_part.period_of_machine;
delta3 = 2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
slot_max = parameters_of_stator.number_of_slot/q;
N_total = 4*(size_of_space_harmonics)+2*slot_max*(lamda_max+1);

%% Generate the matrixPN
% Initial some parameters
matrixPN = zeros(N_total, N_total);
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

%% Boundary conditions in the matrixPN

%1. for boundary condition r = r1;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset1 = 0;
    m = m_offset1+n;
    matrixPN(m, :) = [zeros(1,(n-1)*2), -(space_harmonics*q)*r1^(-(space_harmonics*q)-1), (space_harmonics*q)*r1^((space_harmonics*q)-1), zeros(1, N_total-2*n)];
end

%2. for boundary condition1 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset2 = size_of_space_harmonics;
    m = m_offset2+n;
    matrixPN(m, :) = [zeros(1,(n-1)*2), -(space_harmonics*q)*r2^(-(space_harmonics*q)-1), (space_harmonics*q)*r2^((space_harmonics*q)-1), zeros(1, 2*(size_of_space_harmonics)-2*n),...
        zeros(1,(n-1)*2), mur*(space_harmonics*q)*r2^(-(space_harmonics*q)-1), -mur*(space_harmonics*q)*r2^(space_harmonics*q-1), zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
end

%3. for boundary condition2 r = r2 ;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset3 = 2*(size_of_space_harmonics);
    m = m_offset3+n;
    matrixPN(m, :) = [zeros(1,(n-1)*2), r2^(-(space_harmonics*q)), r2^(space_harmonics*q), zeros(1, 2*(size_of_space_harmonics)-2*n),...
        zeros(1,(n-1)*2), -r2^(-(space_harmonics*q)), -r2^(space_harmonics*q), zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
end


%4. for the boudary condition r = r4;
for slot = slot_min: dslot: slot_max
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset4 = 3*(size_of_space_harmonics);
        m = m_offset4+(slot-1)*(lamda_max+1)+n;
        if(lamda ~= 0)
            matrixPN(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -(lamda*pi/delta3)*r4^(-(lamda*pi/delta3)-1), (lamda*pi/delta3)*r4^((lamda*pi/delta3)-1), zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
        elseif(lamda == 0)
            matrixPN(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                0, 1/r4, zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
        end
    end
end

%5. for the boudary condition1 r = r3;
for slot = slot_min: dslot: slot_max
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset5 = 3*(size_of_space_harmonics)+slot_max*(lamda_max+1);
        m = m_offset5+(slot-1)*(lamda_max+1)+n;
        if(lamda ~= 0)
            matrixPN(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -r3^(-(lamda*pi/delta3)), -r3^(lamda*pi/delta3), zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];    
            h = 0;
            for space_harmonics = space_harmonics_vector
                h = h+1;
                matrixPN(m, :) = matrixPN(m, :)+[zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(h-1)),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*r3^(-space_harmonics*q), f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*r3^(space_harmonics*q),...
                    zeros(1, N_total-2*(size_of_space_harmonics)-2*h)];
            end
        elseif(lamda == 0)
            matrixPN(m, :) = [zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(n-1)), ...
                -1, -log(r3), zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*n)];
            h = 0;
            for space_harmonics = space_harmonics_vector
                h = h+1;
                matrixPN(m, :) = matrixPN(m, :)+[zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(h-1)),...
                    f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*r3^(-space_harmonics*q), f2(lamda)*(M2(space_harmonics, lamda, slot)+M2(space_harmonics, -lamda, slot))*r3^(space_harmonics*q),...
                    zeros(1, N_total-2*(size_of_space_harmonics)-2*h)];
            end
        end
    end
end

% 6. for the boudary condition2 r = r3;
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset6 = 3*(size_of_space_harmonics)+2*slot_max*(lamda_max+1);
    m = m_offset6+n;
    matrixPN(m, :) = [zeros(1, 2*(size_of_space_harmonics)), zeros(1, 2*(n-1)),...
        (-space_harmonics*q)*r3^(-space_harmonics*q-1), (space_harmonics*q)*r3^(space_harmonics*q-1),...
        zeros(1, N_total-2*(size_of_space_harmonics)-2*n)];
    for slot = slot_min: dslot: slot_max
        for lamda = lamda_min: dlamda: lamda_max
            l = round((lamda-lamda_min)/dlamda+1);
            if(lamda ~= 0)
                matrixPN(m, :) = matrixPN(m, :)+[zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(l-1)),...
                    -q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*(-lamda*pi/delta3)*r3^(-lamda*pi/delta3-1),...
                    -q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*(lamda*pi/delta3)*r3^(lamda*pi/delta3-1),...
                    zeros(1, N_total-4*(size_of_space_harmonics)-2*(slot-1)*(lamda_max+1)-2*l)];
            elseif(lamda == 0)
                matrixPN(m, :) = matrixPN(m, :)+[zeros(1, 4*(size_of_space_harmonics)), zeros(1, 2*(slot-1)*(lamda_max+1)), zeros(1, 2*(l-1)),...
                    0, -q/(4*pi)*(M1(space_harmonics, lamda, slot)+M1(space_harmonics, -lamda, slot))*1/r3,...
                    zeros(1, N_total-4*(size_of_space_harmonics)- 2*(slot-1)*(lamda_max+1)-2*l)];
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

