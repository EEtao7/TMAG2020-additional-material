function [B_radial, B_tangential] = Magnetic_Field_Noload_Slotless_time_stepping(parameters_of_stator, parameters_of_rotor, parameters_of_space_harmonics, parameters_of_other_part, thetam, thetam_shift)
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
harmonics_max = parameters_of_space_harmonics.PM_space_harmonics_max;

%% Calculate the magnetic field
thetam = thetam-thetam_shift; % rotor rotates with mechanical degree shift (time)
B_radial  = 0;
B_tangential = 0;
if(r<Rs && r>Rm)
    for i =1: 1: harmonics_max
        n = 2*i-1;
        % when np ~= 1
        B_radial = B_radial + (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
            ((r/Rs)^(n*Pr-1)*(Rm/Rs)^(n*Pr+1)+(Rm/r)^(n*Pr+1))...
            *(n*Pr-1+2*(Rr/Rm)^(n*Pr+1)-(n*Pr+1)*(Rr/Rm)^(2*n*Pr))/((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))*cos(n*Pr*thetam);
        B_tangential = B_tangential + (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
            (-(r/Rs)^(n*Pr-1)*(Rm/Rs)^(n*Pr+1)+(Rm/r)^(n*Pr+1))...
            *(n*Pr-1+2*(Rr/Rm)^(n*Pr+1)-(n*Pr+1)*(Rr/Rm)^(2*n*Pr))/((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))*sin(n*Pr*thetam);
    end
elseif(r>Rr && r<=Rm)
    for i =1: 1: harmonics_max
        n = 2*i-1;
        % when np ~= 1
        B_radial = B_radial + (Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
            ((Rr/Rm)^(n*Pr-1)*(Rr/r)^(n*Pr+1)+(r/Rm)^(n*Pr-1))...
            *((n*Pr-1/mur)*(Rm/Rs)^(2*n*Pr)+(1+1/mur)*(Rr/Rm)^(n*Pr+1)*(Rm/Rs)^(2*n*Pr)-(n*Pr+1/mur)-(1-1/mur)*(Rr/Rm)^(n*Pr+1))/...
            ((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))*cos(n*Pr*thetam)+...
            ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1)))*(Rr/r)^(n*Pr+1)*cos(n*Pr*thetam)+...
            ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*((n*Pr)^2/((n*Pr)^2-1)))*cos(n*Pr*thetam);
        B_tangential = B_tangential + (-Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1))*...
            (-(Rr/Rm)^(n*Pr-1)*(Rr/Rm)^(n*Pr+1)+(r/Rm)^(n*Pr-1))...
            *((n*Pr-1/mur)*(Rm/Rs)^(2*n*Pr)+(1+1/mur)*(Rr/r)^(n*Pr+1)*(Rm/Rs)^(2*n*Pr)-(n*Pr+1/mur)-(1-1/mur)*(Rr/Rm)^(n*Pr+1))/...
            ((mur+1)/mur*(1-(Rr/Rs)^(2*n*Pr))-(mur-1)/mur*((Rm/Rs)^(2*n*Pr)-(Rr/Rm)^(2*n*Pr)))*sin(n*Pr*thetam)+...
            ((Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*(n*Pr/((n*Pr)^2-1)))*(Rr/r)^(n*Pr+1)*sin(n*Pr*thetam)+...
            ((-Br/mur)*(4/n/pi)*sin(n*pi*alpha/2)*((n*Pr)^1/((n*Pr)^2-1)))*sin(n*Pr*thetam);
    end
else
    disp('Beyond the magnetic field solution domain.');
end
end

