function [B_radialn, B_tangentialn] = Amplitude_of_B_PM_time_harmonics(parameters_of_stator, parameters_of_rotor,  parameters_of_other_part, parameters_of_space_harmonics, time_harmonics)
%Magnetic_Field_Noload_Slotless is used to calculate the magnetic field in
% the slotless PMSM without the current (just consider the PM)
%% Define some useful parameters
Br = parameters_of_rotor.B_remanence;
mur = parameters_of_rotor.mur;
alpha = parameters_of_rotor.ratio_of_magnet;
Pr = parameters_of_rotor.pole_pairs_of_rotor;
Rs = parameters_of_stator.radius_of_stator;
Rr = parameters_of_rotor.radius_of_rotor;
Rm = parameters_of_rotor.radius_of_rotor+parameters_of_rotor.thickness_of_PM;
r = parameters_of_other_part.radius;
space_harmonics_vector = parameters_of_space_harmonics.PM_space_harmonics;

%% Calculate the magnetic field
size_of_space_harmonics = length(space_harmonics_vector);
B_radialn  = zeros(1, size_of_space_harmonics);
B_tangentialn = zeros(1, size_of_space_harmonics);

if(r<Rs && r>Rm)
    i = 0;
    for space_harmonics =space_harmonics_vector
        i = i+1;
        n = space_harmonics/Pr;
        if(space_harmonics == time_harmonics)
            % when np ~= 1
            B_radialn(i) =  (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
                ((r/Rs)^(n*Pr-1)*(Rm/Rs)^(n*Pr+1)+(Rm/r)^(n*Pr+1))...
                *(n*Pr-1+2*(Rr/Rm)^(n*Pr+1)-(n*Pr+1)*(Rr/Rm)^(2*n*Pr))/((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)));
            B_tangentialn(i) = (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
                (-(r/Rs)^(n*Pr-1)*(Rm/Rs)^(n*Pr+1)+(Rm/r)^(n*Pr+1))...
                *(n*Pr-1+2*(Rr/Rm)^(n*Pr+1)-(n*Pr+1)*(Rr/Rm)^(2*n*Pr))/((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)));
        else
            B_radialn(i) = 0;
            B_tangentialn(i) = 0;
        end
    end
elseif(r>Rr && r<=Rm)
    i = 0;
    for space_harmonics =space_harmonics_vector
        i = i+1;
        n = space_harmonics/Pr;
        if(space_harmonics == time_harmonics)
            % when np ~= 1
            B_radialn(i) = (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
                ((Rr/Rm)^(n*Pr-1)*(Rr/r)^(n*Pr+1)+(r/Rm)^(n*Pr-1))...
                *((n*Pr-1/mur)*(Rm/Rs)^(2*n*Pr)+(1+1/mur)*(Rr/Rm)^(n*Pr+1)*(Rm/Rs)^(2*n*Pr)-(n*Pr+1/mur)-(1-1/mur)*(Rr/Rm)^(n*Pr+1))/...
                ((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))+...
                ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1)))*(Rr/r)^(n*Pr+1)+...
                ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*((n*Pr)^2/((n*Pr)^2-1)));
            B_tangentialn(i) = (-Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
                (-(Rr/Rm)^(n*Pr-1)*(Rr/Rm)^(n*Pr+1)+(r/Rm)^(n*Pr-1))...
                *((n*Pr-1/mur)*(Rm/Rs)^(2*n*Pr)+(1+1/mur)*(Rr/r)^(n*Pr+1)*(Rm/Rs)^(2*n*Pr)-(n*Pr+1/mur)-(1-1/mur)*(Rr/Rm)^(n*Pr+1))/...
                ((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))+...
                ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1)))*(Rr/r)^(n*Pr+1)+...
                ((-Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*((n*Pr)^1/((n*Pr)^2-1)));
        else
            B_radialn(i) = 0;
            B_tangentialn(i) = 0;
        end
    end
else
    disp('Beyond the magnetic field solution domain.');
end
end

