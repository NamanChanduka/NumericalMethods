%
% -------------------
% This is a solution to the Homework no. 5 for  CE 465, Numerical Methods in Civil Engineering 
% course at IIT Bombay instructed by Prof. Ravi Sinha (Spring 2020-21)
% 
% in the present file we have to evaluate the roots of the functions
% f1(x) = x^2 + y^2 - 4 = 0 and 
% f2(x) = x^2 - y^2 - 1 = 0 using newtons method
% 
% 
% Author: Naman Chanduka, IIT Bombay
%
%
% -------------------

clear; clc; close all; %used for clearing the workspace

tic; %used for monitoring the runtime of the program

format long; %to accurately see more decimal points

%%
%%%%%%%%  Start of User Input  %%%%%%%%%%%%%%%
 
x0 = 1.6;                       %initial guess for x

y0 = 1.2;                       %initial guess for y

conv = 2.22e-14;                %convergence criterion

%%%%%%%% End of User Input  %%%%%%%%%%%%%


%%
%%%% Start of Program %%%%
F1 = @(x) x(1).^2 + x(2).^2 - 4;                                                                %Function 1 defined
F2 = @(x) x(1).^2 - x(2).^2 - 1;                                                                %Function 2 defined

x_now = [x0;y0];                                                                                %initialize x_now
x_next = [0;0];
iteration = 0;                                                                                  %initialize iterations

looping_condition = abs(x_next(1) - x_now(1)) > conv || abs(x_next(2) - x_now(2)) > conv;       %looping condition created to start while loop
while (looping_condition)
    if(x_now(1) == x0 && x_next(1) == 0 && x_now(2) == y0 && x_next(2) == 0)                    %if first iteration then dont update x_now, use initialized value
        
    else
        x_now = x_next;                                                                         %if not first iteration, then update x_now
    end
    Jacobian = [2*x_now(1), 2*x_now(2); 2*x_now(1), -2*x_now(2)];                               %Jacobian [2x, 2y; 2x, -2y] 
    F = [F1(x_now); F2(x_now)];                                                                 %function vector [f1;f2]
    delta = [ sum(F)/(-4*x_now(1)) ; (F2(x_now)-F1(x_now))/(4*x_now(2)) ] ; 
    x_next = x_now + delta;                                                                     %Delta vector is calculated so that there is no need to calculate inverse of Jacobian, so higher accuracy is preserved
    iteration = iteration + 1;                                                                  %updated iteration variable
    y1(iteration) = abs(x_next(1)-x_now(1));
    y2(iteration) = abs(x_next(2)-x_now(2));
    x(iteration) = iteration;
    looping_condition = abs(x_next(1) - x_now(1))/x_now(1) > conv || abs(x_next(2) - x_now(2))/x_now(2) > conv;   %updated looping condition to continue while loop
 
end

figure(1);    
semilogy(x, y1);                                                               %plot for the relative error root x.
grid on;
xlabel('Number of iterations');
ylabel('Relative error (log scale)');
title('Newton method');

figure(2);    
semilogy(x, y2);                                                               %plot for the relative error root y.
grid on;
xlabel('Number of iterations');
ylabel('Relative error (log scale)');
title('Newton method');

disp(x_next);
disp(iteration);
toc;
%%%% End of Program %%%%