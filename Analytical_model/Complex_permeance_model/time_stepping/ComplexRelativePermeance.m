function lamda= ComplexRelativePermeance(parameters_of_stator, parameters_of_rotor, complex_coordinate)
%ComplexRelativePermeance is used to generate the complex relative
%permeance in the corresponding complex coordinate;
%% Define some useful parameters
Rs = parameters_of_stator.radius_of_stator;
Rr = parameters_of_rotor.radius_of_rotor;
Ns = parameters_of_stator.number_of_slot;
alpha = parameters_of_stator.ratio_of_tooth;
s = complex_coordinate;
deg = pi/180;

%% Calculate some varible in the calculation process
thetas = (360/Ns)*deg;
theta2 = (180/Ns*(1-alpha)+180/Ns)*deg;
b0 = ((1-alpha)*360/Ns)*deg;
g = log(Rs/Rr);
b = (b0/(2*g)+sqrt((b0/(2*g))^2+1))^2;
a = 1/b;

%% Solve the w
% use the optimization tool in the MATLAB to solve the equations
    function F = fun(y)
        p = sqrt((y-b)/(y-a));
        F = 1i*g/pi*(log((1+p)*(b-p)/(1-p)/(b+p))-2*(b-1)/sqrt(b)*atan(p/sqrt(b)))+log(Rs)+1i*theta2-log(s);
    end
x0 = 0+1i;
[w, fval] = fsolve(@fun, x0);

% If the result has big error, the program will use the VPAsolve in lower
% speed
if(abs(fval)>1e-5)
    w = VPAsolve(a, b, g, Rs, theta2, s);
end

% check the precision of result
% p = sqrt((w-b)/(w-a));
% result = 1i*g/pi*(log((1+p)*(b-p)/(1-p)/(b+p))-2*(b-1)/sqrt(b)*atan(p/sqrt(b)))+log(Rs)+1i*theta2-log(s);

% syms x
% p = sqrt((x-b)/(x-a));
% f = 1i*g/pi*(log((1+p)*(b-p)/(1-p)/(b+p))-2*(b-1)/sqrt(b)*atan(p/sqrt(b)))+log(Rs)+1i*theta2-log(s);
% w = vpasolve(f, x);

%% Calculate the lamda
k = Rs*exp(1i*(g/pi*log(w)+thetas*0.5));
lamda = k/s*(w-1)/sqrt(w-a)/sqrt(w-b);

end

