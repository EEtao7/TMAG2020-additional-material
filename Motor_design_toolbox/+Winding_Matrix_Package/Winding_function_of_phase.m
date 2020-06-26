function N_of_phase = Winding_function_of_phase(parameters_of_stator, thetam, windingmatrix_of_phase)
% Winding_function_of_phase is used to generate the winding function of one phase based on the winding matrix of one phase
%% Define some parameters
Ns = parameters_of_stator.number_of_slots;
Np = parameters_of_stator.number_of_phase;
turns = parameters_of_stator.turns_of_phase*2/(Ns*2/Np); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数
y = parameters_of_stator.pitch_of_coils;
thetam_min = thetam.thetam_min;
dthetam = thetam.dthetam;
thetam_max = thetam.thetam_max;
size_of_theta = thetam.size_of_thetam;

number_of_slot_per_phase = Ns/Np;
alpha = 360/Ns;

%% Generate Winding function of phase
% Check the Winding Matrix of Phase
[~, size_of_matrix] = size(windingmatrix_of_phase);
if(number_of_slot_per_phase == round(size_of_matrix/2))
    %% Generate the winding function of every coils
    N_assistant_matrix = zeros(number_of_slot_per_phase, size_of_theta);
    for m=1 : 1: number_of_slot_per_phase
        theta1 = (abs(windingmatrix_of_phase(2*m-1))-1)*alpha;
        theta2 = (abs(windingmatrix_of_phase(2*m))-1)*alpha;
        if(theta1<theta2)
            for theta = thetam_min: dthetam: thetam_max
                n = round((theta-thetam_min)/dthetam+1);
                if(theta>=theta1 && theta<theta2)
                    N_assistant_matrix(m, n) = turns*(1-y/Ns)*sign(windingmatrix_of_phase(2*m-1));          
                else
                    N_assistant_matrix(m, n) = -turns*(y/Ns)*sign(windingmatrix_of_phase(2*m-1));   
                end
            end
        elseif(theta1>theta2)
            for theta = thetam_min: dthetam: thetam_max
                n = round((theta-thetam_min)/dthetam+1);
                if(theta>=theta1 || theta<theta2)
                    N_assistant_matrix(m, n) = turns*(1-y/Ns)*sign(windingmatrix_of_phase(2*m-1));          
                else
                    N_assistant_matrix(m, n) = -turns*(y/Ns)*sign(windingmatrix_of_phase(2*m-1));   
                end
            end
        else
            disp('elements errors of Winding Matrix of phase');
        end
    end
    %% Combine the winding function of every coils into a winding function of one phase
    N_of_phase = sum(N_assistant_matrix);
else
    disp('size errors of Winding Matrix of phase')
     N_of_phase = zeros(1, size_of_theta);
end

end

