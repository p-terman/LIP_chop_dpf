function [a, b, sigma_a, sigma_b xi2_v] = linfit(X,Y,error)
% function [a, b, sigma_a, sigma_b xi2_v] = linfit(X,Y,error)
%
% This function accepts vectors X, Y and an error for Y, and returns the
% coefficients a and b for the equation Y = a + bX with respective errors
% and reduced chi squared
%
% CHF

N = length(X);
S = 0; Sx = 0; Sy = 0; Stt = 0; b = 0;

for i = 1:N
    S = S + (1/error(i))^2;
    Sx = Sx + (X(i)/error(i)^2);
    Sy = Sy + (Y(i)/error(i)^2);
end

for i = 1:N
    t(i) = (X(i) - (Sx/S))/error(i);
    Stt = Stt + t(i)^2;
    b = b + (t(i)*Y(i)/error(i));
end

b = b/Stt;
a = (Sy - Sx*b)/S;

sigma_a = sqrt((1 + (Sx^2/S/Stt))/S);
sigma_b = sqrt(1/Stt);
Covab = -Sx/S/Stt;
rab = Covab/sigma_a/sigma_b;


if nargout > 4
    s = (Y - (a + b.*X))./error; %difference over sigma
    xi2 = sum(s.^2); %compute xi^2
    v = numel(Y) - 2; %number of degrees of freedom
    xi2_v = xi2/v; %reduce
end