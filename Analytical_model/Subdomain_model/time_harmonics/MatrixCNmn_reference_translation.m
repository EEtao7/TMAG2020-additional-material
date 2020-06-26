function matrixCNmn_single = MatrixCNmn_reference_translation(parameters_of_stator, parameters_of_other_part, parameters_of_space_harmonics, matrixCNmn_single_ref)
% This function is used to transfer the matrixCNmn_based on the radious
% reference to matrixCNmn_single in the fact

%% Define some useful parameters
unit_coefficient = parameters_of_other_part.unit_coefficient;
r_ref1 = parameters_of_other_part.r_reference1/unit_coefficient;
r_ref2 = parameters_of_other_part.r_reference2/unit_coefficient;
q = parameters_of_other_part.period_of_machine;
delta3 = 2*pi*(1-parameters_of_stator.ratio_of_tooth)/parameters_of_stator.number_of_slot;

space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection/q;
size_of_space_harmonics = length(space_harmonics_vector);
lamda_max = parameters_of_space_harmonics.lamda_max;
unit_slot = parameters_of_stator.unit_slots;
N_total = 4*(size_of_space_harmonics)+2*unit_slot*(lamda_max+1);

%% Generate the translation matrix T
% Initial some parameters
matrix_T = zeros(N_total, 1);
lamda_min = 0;
dlamda = 1;
slot_min = 1;
dslot = 1;

% for the C11, C12
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset1 = 0;
    m = m_offset1+n;
    matrix_T(round(2*m-1), 1) = r_ref1^(space_harmonics*q);
    matrix_T(round(2*m), 1) = 1/r_ref1^(space_harmonics*q);
end

% for the C21, C22
n = 0;
for space_harmonics = space_harmonics_vector
    n = n+1;
    m_offset2 = size_of_space_harmonics;
    m = m_offset2+n;
    matrix_T(round(2*m-1), 1) = r_ref1^(space_harmonics*q);
    matrix_T(round(2*m), 1) = 1/r_ref1^(space_harmonics*q);
end

%4. for the C31(s, lamda), C32(s, lamda)
for slot = slot_min: dslot: unit_slot
    for lamda = lamda_min: dlamda: lamda_max
        n = round((lamda-lamda_min)/dlamda+1);
        m_offset4 = 2*(size_of_space_harmonics);
        m = m_offset4+(slot-1)*(lamda_max+1)+n;
        matrix_T(round(2*m-1), 1) = r_ref2^(lamda*pi/delta3);
        matrix_T(round(2*m), 1) = 1/r_ref2^(lamda*pi/delta3);
    end
end

%% Calculate the matrixCNmn_single
matrixCNmn_single = matrixCNmn_single_ref.*matrix_T;

end

