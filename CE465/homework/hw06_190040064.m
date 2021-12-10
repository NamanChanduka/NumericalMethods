%
% -------------------
% This is a solution to the Homework no. 6 for  CE 465, Numerical Methods in Civil Engineering
% course at IIT Bombay instructed by Prof. Ravi Sinha (Spring 2020-21)
%
% in the present file we have to evaluate the eight storied industrial
% structure as given in the question
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

k = [100; 100; 100; 80; 80; 1; 60; 60; 20; 20; 0.15; 0.2];                  %Defining the stiffness values

A = [k(1)+k(2)+k(10), -k(2), -k(10), 0, 0, 0, 0, 0;                         %Defining stiffness matrix which is A matrix in Ax = b in this case
    -k(2), k(2)+k(3)+k(11), -k(3), 0, 0, 0, 0, -k(11);
    -k(10), -k(3), k(3)+k(4)+k(9)+k(10), -k(4), -k(9), 0, 0, 0;
    0, 0, -k(4), k(4)+k(5)+k(12), -k(5), 0, -k(12), 0;
    0, 0, -k(9), -k(5), k(5)+k(6)+k(9), -k(6), 0, 0;
    0, 0, 0, 0, -k(6), k(6)+k(7), -k(7), 0;
    0, 0, 0, -k(12), 0, -k(7), k(7)+k(8)+k(12), -k(8);
    0, -k(11), 0, 0, 0, 0, -k(8), k(8)+k(11)];

b = [0; 0; 0; 30; 0; 20; 0; 20];                                            %Defining force vector which is b vector in Ax = b in this case

n = 8;                                                                      %initialise number of equations according to original matrices provided

%%%%%%%% End of User Input  %%%%%%%%%%%%%


%%
%%%% Start of Program %%%%

pivot = double.empty(n, 0);                                                 %initialise pivot as 0 vector initially

[A, pivot, determinant, ier] = hw06_factor(A, n, pivot);                    %storing new factorised matrix instead of original matrix to reduce space complexity

if ier == 1  %checking if error flag equals 1
    disp('Solution Doesnt exist')
elseif ier == 0
    b = hw06_solve(A, n, b, pivot);                                         %storing answer (x vector) in b vector from Ax = b to reduce space complexity
    disp('Displacement at each corresponding DOF: ');
    disp(b);
elseif ier == -1
    disp("Some error has occured during factorisation");
end

toc;
%%%% End of Program %%%%

