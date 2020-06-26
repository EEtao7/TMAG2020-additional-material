function w = VPAsolve(a, b, g, Rs, theta2, s)
%VPAsolve is used to calculate the root of the function
syms x
p = sqrt((x-b)/(x-a));
f = 1i*g/pi*(log((1+p)*(b-p)/(1-p)/(b+p))-2*(b-1)/sqrt(b)*atan(p/sqrt(b)))+log(Rs)+1i*theta2-log(s);
w = vpasolve(f, x);
end

