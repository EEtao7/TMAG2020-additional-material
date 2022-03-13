function matrixCNmn_negative_time_harmonics = MatrixCNmn_negative_time_harmonics(parameters_of_space_harmonics, matrixCNmn_single)
% This function is used to expand the matrixCNmn
%% Define some useful parameters
space_harmonics_vector = parameters_of_space_harmonics.space_harmonics_collection;
size_of_space_harmonics_vector = length(space_harmonics_vector);
matrixCNmn_positive_time_harmonics = matrixCNmn_single;

%% Expand the matrixCNmn when use negative time harmonics
vector1 = conj(matrixCNmn_positive_time_harmonics( 2*size_of_space_harmonics_vector: -1: size_of_space_harmonics_vector+1 ))';
vector2 = conj(matrixCNmn_positive_time_harmonics( size_of_space_harmonics_vector: -1: 1 ))';
vector3 = conj(matrixCNmn_positive_time_harmonics(4*size_of_space_harmonics_vector: -1: 3*size_of_space_harmonics_vector+1))';
vector4 = conj(matrixCNmn_positive_time_harmonics(3*size_of_space_harmonics_vector: -1: 2*size_of_space_harmonics_vector+1))';
vector5 = conj(matrixCNmn_positive_time_harmonics(4*size_of_space_harmonics_vector+1: end))';
matrixCNmn_negative_time_harmonics = [vector1, vector2, vector3, vector4, vector5];
end

