function Brn = Br_radial_megnetization(parameters_of_rotor, harmonics_q)
% caculate the amplititude of remanence Brn in the radial megenetization 
p = parameters_of_rotor.pole_pairs_of_rotor;
n = harmonics_q;

if(mod((n/p), 2) == 1)
    Brn = parameters_of_rotor.B_remanence*2*p*...
        sin(n*pi*parameters_of_rotor.ratio_of_magnet/(2*p))/(pi*n);
else
    Brn = 0;
end

end