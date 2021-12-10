%
% -------------------
% This is a solution to the Homework no. 7 for  CE 465, Numerical Methods in Civil Engineering
% course at IIT Bombay instructed by Prof. Ravi Sinha (Spring 2020-21)
%
% in the present file we have to evaluate the solution for the given system
% of equations using 2 different methods, Gauss Elim with pivoting and
% Conjugate Gradient Method
%
% Author: Naman Chanduka, IIT Bombay
%
%
% -------------------

clear; clc; close all; %used for clearing the workspace

format long; %to accurately see more decimal points

%%
%%%%%%%%  Start of User Input  %%%%%%%%%%%%%%%


n = 5;                                                                      %initialise size of hilbert matrix and number of equations

kmax = 5;                                                                   %kmax for conjugate gradient method

%%%%%%%% End of User Input  %%%%%%%%%%%%%


%%
%%%% Start of Program %%%%

for i = 1:n                                                                 %Defining b, H and H_inverse
    b(i) = (-1)^(i+1);
    for j = 1:n
        H(i,j) = 1/(i+j-1);
        H_inverse(i,j) = (((-1)^(i+j))*factorial(n+i-1)*factorial(n+j-1))/((i+j-1)*((factorial(i-1)*factorial(j-1))^2)*factorial(n-i)*factorial(n-j));
    end
end

actual_answer = H_inverse*(b');                                             %Actual answer for comparison using the inverse of Hilbert matrix

pivot = double.empty(n, 0);                                                 %initialise pivot as 0 vector initially

timeForGauss = tic;

[A, pivot, determinant, ier] = hw07_factor(H, n, pivot);                    %storing new factorised matrix

if ier == 1  %checking if error flag equals 1
    disp('Solution Doesnt exist')
elseif ier == 0
    x_gauss = hw07_solve(A, n, b, pivot);                                   %storing answer (x vector) 
    disp(x_gauss');
elseif ier == -1
    disp("Some error has occured during factorisation");
end

toc(timeForGauss);

timeForCGM = tic;

x_CGM = hw07_conjugate(H, b', n, 100*eps, kmax);
disp(x_CGM);

toc(timeForCGM);

disp(actual_answer);

%%%% End of Program %%%%


