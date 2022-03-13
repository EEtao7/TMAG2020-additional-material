function n2 = Current_space_harmonics(parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, parameters_of_time_harmonics)
% This function is used to generate the space harmonics with the Current source
%% Define some useful parameters
c2_max = parameters_of_space_harmonics.space_harmonics_constant_max;
hc = parameters_of_time_harmonics.current_time_harmonics;
Ns = parameters_of_stator.number_of_slot;
t = parameters_of_other_part.period_of_machine;
z = Z(Ns,t);
Np = parameters_of_stator.number_of_phase;
n2 = [];

%% Generate the space harmonics with the PM
if(sum(hc)==0)
    disp('motor works in the Noload situation');
else
    for c2 = -c2_max: 1: c2_max
        n_assistant1 = hc+c2*(z*Np*t);
        n_assistant2 = -hc+c2*(z*Np*t);
        n2 = union(n2, union(n_assistant1, n_assistant2));
    end
end

    function z = Z(Ns, t)
        if(mod(round(Ns/t), 2)==1)
            z = 1;
        elseif(mod(round(Ns/t), 2)==0)
            z = 2;
        end
    end
end

