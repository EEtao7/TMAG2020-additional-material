function unit_windingmatrix_of_phaseA= Unit_WindingMatrix_of_PhaseA(parameters_of_windings)
% The function Unit_WindingMatrix_of_PhaseA is used to generate the winding matrix of phase A in the unit machine period
%%  Define some parameters
Ns = parameters_of_windings.number_of_slots;
Ps = parameters_of_windings.pole_pairs_of_stator;
Np = parameters_of_windings.number_of_phase;
turns = parameters_of_windings.turns_of_phase*2/(Ns*2/Np); % 第一个×2表示匝数换算成导体数，第二个×2代表双层绕组的绕组层数

t = gcd(Ns, Ps); % period of machine
y = round(Ns/(2*Ps)); % pitch of coils classical setting
number_of_spokes = Ns/t;
angle_of_spoke = 360/number_of_spokes;
number_of_spokes_per_phase = number_of_spokes/Np;
relative_tolerance = 1e-6*360;

%% Generate the Unit winding matrix of PhaseA
% Check the slot and pole-pairs combination
if (mod(Ns,Np*t)==0)
    disp('The Windings satify the slots/poles combination');
    
    %% Generate the Matrix of eletrical degree and sequence
    dmechanical_angle = angle_of_spoke;
    mechanical_angle_min = 0;
    mechanical_angle_max = 360-dmechanical_angle;
    mechanical_angle = mechanical_angle_min : dmechanical_angle : mechanical_angle_max;
    electrical_angle = mod(mechanical_angle*(Ps/t), 360);
    Matrix_A = zeros(2,number_of_spokes); %  Matrix of eletrical degree and sequence
    Matrix_A(1, :) = 1: 1: number_of_spokes;
    Matrix_A(2, :) = electrical_angle;
    
    %% Generate the Winding Matrix of Phase A in the unit machine
    % Generate the assistant matrix for the windings
    assistant_matrix = zeros(2, number_of_spokes_per_phase);
    if (mod(number_of_spokes_per_phase, 2)==0)
        for m = 1:1:round(number_of_spokes_per_phase/2)
            col1 = find(abs(Matrix_A(2, :) - mod(angle_of_spoke*(m-1),360))<relative_tolerance);
            col2 = find(abs(Matrix_A(2, :) - mod(180+angle_of_spoke*(m-1),360))<relative_tolerance);
            assistant_matrix(1, m) = Matrix_A(1, col1);
            assistant_matrix(2, m) = 1;
            if(m+round(number_of_spokes_per_phase/2)<=number_of_spokes_per_phase)
                assistant_matrix(1, m+round(number_of_spokes_per_phase/2)) = Matrix_A(1, col2);
                assistant_matrix(2, m+round(number_of_spokes_per_phase/2)) = -1;
            end
        end
    elseif(mod(number_of_spokes_per_phase, 2)==1)
        for m = 1:1:ceil(number_of_spokes_per_phase/2)
            col1 = find(abs(Matrix_A(2, :) - mod(angle_of_spoke*(m-1),360))<relative_tolerance);
            col2 = find(abs(Matrix_A(2, :) - mod(180+angle_of_spoke*(m-1)+angle_of_spoke/2,360))<relative_tolerance);
            assistant_matrix(1, m) = Matrix_A(1, col1);
            assistant_matrix(2, m) = 1;
            if(m+ceil(number_of_spokes_per_phase/2)<=number_of_spokes_per_phase)
                assistant_matrix(1, m+ceil(number_of_spokes_per_phase/2)) = Matrix_A(1, col2);
                assistant_matrix(2, m+ceil(number_of_spokes_per_phase/2)) = -1;
            end
        end
    end
    
    % Generate the Winding matrix of Phase A in the unit machine
    unit_windingmatrix_of_phaseA = zeros(3, 2*number_of_spokes_per_phase);
    for m = 1 : 1 : number_of_spokes_per_phase
        unit_windingmatrix_of_phaseA(1, 2*m-1) = assistant_matrix(1,m);
        unit_windingmatrix_of_phaseA(2, 2*m-1) = 2;
        unit_windingmatrix_of_phaseA(3, 2*m-1) = assistant_matrix(2,m)*turns;
        unit_windingmatrix_of_phaseA(1, 2*m) = mod(assistant_matrix(1,m)+y, Ns);
        unit_windingmatrix_of_phaseA(2, 2*m) = 1;
        unit_windingmatrix_of_phaseA(3, 2*m) = -1*assistant_matrix(2,m)*turns;
    end
    unit_windingmatrix_of_phaseA(1, :) = round(unit_windingmatrix_of_phaseA(1, :));
    
    % Do some correction of some elements in the winding Matrix
    % corrent 0 elements in the winding matrix
    [row, col] = find(unit_windingmatrix_of_phaseA == 0);
    [size_of_row, ~] = size(row);
    for m = 1: size_of_row
        unit_windingmatrix_of_phaseA(row(m), col(m)) = abs((unit_windingmatrix_of_phaseA(row(m), col(m))+Ns));
    end
    
else
    disp('slots and pole-pairs combination have errors');
    unit_windingmatrix_of_phaseA = zeros(1,1);
end


end

