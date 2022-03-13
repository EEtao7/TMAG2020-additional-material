function Brmn = Br_radial_megnetization_time_harmonics(parameters_of_rotor, parameters_of_time_harmonics, harmonics_q, time_harmonics)
% caculate the amplititude of remanence Brn in the radial megenetization 
pr = parameters_of_rotor.pole_pairs_of_rotor;
hm = union(-parameters_of_time_harmonics.PM_time_harmonics, parameters_of_time_harmonics.PM_time_harmonics);
% hm =  parameters_of_time_harmonics.PM_time_harmonics;
m = time_harmonics;
n = harmonics_q;
Brmn = 0;

if(ismember(m, hm))% 判断时间谐波是否为以PM为激励源的时间谐波
    if(round(n)==round(m))% 判断空间谐波是否等于时间谐波
        if(mod((n/pr), 2) == 1)% 判断空间谐波是否为奇数
            Brmn = parameters_of_rotor.B_remanence*2*pr*...
                sin(n*pi*parameters_of_rotor.ratio_of_magnet/(2*pr))/(pi*n);
        end
    end
end
end