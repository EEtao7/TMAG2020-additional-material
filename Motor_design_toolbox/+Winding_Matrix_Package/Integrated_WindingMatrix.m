function [integrated_windingmatrix_of_phaseA, integrated_windingmatrix] =Integrated_WindingMatrix (parameters_of_stator, Unit_WindingMatrix_of_PhaseA)
%The function Integrated_WindingMatrix is used to generate the Integrated Winding Matrix of phase A and Winding Matrix for the whole machine
%% Define some parameters
Ns = parameters_of_stator.number_of_slots;
Ps = parameters_of_stator.pole_pairs_of_stator;
Np = parameters_of_stator.number_of_phase;

t = gcd(Ns, Ps); % period of machine
number_of_spokes_per_phase = Ns/Np/t;

%% Expand the Unit A Phase Matrix to the Integrated WindingMatrix of phaseA
integrated_windingmatrix_of_phaseA = zeros(1, 2*Ns/Np);
for m = 1: 1: t
    integrated_windingmatrix_of_phaseA(1, (1: 2*number_of_spokes_per_phase)+(m-1)*2*number_of_spokes_per_phase) = (Unit_WindingMatrix_of_PhaseA(1,:)+(m-1)*round(Ns/t)).*sign(Unit_WindingMatrix_of_PhaseA(3, :));
end

%% Expand the Integrated WindingMatrix of phaseA to Integrated WindingMatrix
integrated_windingmatrix = zeros(Np, 2*Ns/Np);
for n = 1:1:Np
    integrated_windingmatrix(n, :) = mod(abs(integrated_windingmatrix_of_phaseA)+(1/t)*((n-1)/(Np))*Ns, Ns).*sign(integrated_windingmatrix_of_phaseA);
end
% do some correction of some elements in the winding Matrix
[row, col] = find(integrated_windingmatrix == 0);
[size_of_row, ~] = size(row);
for m = 1: size_of_row
    integrated_windingmatrix(row(m), col(m)) = (integrated_windingmatrix(row(m), col(m))+Ns).*sign(integrated_windingmatrix_of_phaseA(1,col(m)));
end
% correct some integrated_windingmatrix_of_phaseA
integrated_windingmatrix_of_phaseA = integrated_windingmatrix(1, :);

end

