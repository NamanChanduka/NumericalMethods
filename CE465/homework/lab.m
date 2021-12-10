format short

clear; clc;
m = 21.5;
h = 11;
d = 60 - 0.7*m;
y0 = sqrt(h^2 + d^2) - d;
for x = 22:1:54 
    x0 = 60-x;
    y(x-21) = sqrt(2*x0*y0 + y0^2);
end
disp(y');