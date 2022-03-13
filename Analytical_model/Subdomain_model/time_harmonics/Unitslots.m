function s = Unitslots(parameters_of_stator, parameters_of_other_part, parameters_of_time_harmonics)
% This function is used to generate the space harmonics with the Current source
%% Define some useful parameters
hc = parameters_of_time_harmonics.current_time_harmonics;
Ns = parameters_of_stator.number_of_slot;
t = parameters_of_other_part.period_of_machine;
z = Z(Ns,t);
Np = parameters_of_stator.number_of_phase;

%% Generate the space harmonics with the PM
if(sum(hc)==0)
    s = 1;
else
    s =  round(Ns/(z*Np*t));
end

    function z = Z(Ns, t)
        if(mod(round(Ns/t), 2)==1)
            z = 1;
        elseif(mod(round(Ns/t), 2)==0)
            z = 2;
        end
    end
end

