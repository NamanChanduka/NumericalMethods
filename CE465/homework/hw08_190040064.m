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

k = [100; 100; 100; 80; 80; 1; 60; 60; 20; 20; 0.15; 0.2];  

K = [k(1)+k(2)+k(10), -k(2), -k(10), 0, 0, 0, 0, 0;                         %Defining stiffness matrix which is A matrix in Ax = b in this case
    -k(2), k(2)+k(3)+k(11), -k(3), 0, 0, 0, 0, -k(11);
    -k(10), -k(3), k(3)+k(4)+k(9)+k(10), -k(4), -k(9), 0, 0, 0;
    0, 0, -k(4), k(4)+k(5)+k(12), -k(5), 0, -k(12), 0;
    0, 0, -k(9), -k(5), k(5)+k(6)+k(9), -k(6), 0, 0;
    0, 0, 0, 0, -k(6), k(6)+k(7), -k(7), 0;
    0, 0, 0, -k(12), 0, -k(7), k(7)+k(8)+k(12), -k(8);
    0, -k(11), 0, 0, 0, 0, -k(8), k(8)+k(11)];

m = [1000, 1000, 1000, 800, 800, 600, 600, 600];                            %Defining the values of m1....m8

n = 8;


%%
%%%% Start of Program %%%%

M_halfinverse = diag([1/sqrt(m(1)), 1/sqrt(m(2)), 1/sqrt(m(3)), 1/sqrt(m(4)), 1/sqrt(m(5)), 1/sqrt(m(6)), 1/sqrt(m(7)), 1/sqrt(m(8))]);         %Creating M^-1/2 matrix

A = M_halfinverse*K*M_halfinverse;                                          %Defining A matrix as specified in the problem statement              

tol = eps;
[eigenval, iterations, offList] = hw08_Jacobi(A, n, tol);                   %Calling Jacobi Function 
iteration = 1:iterations+1;

plot(iteration, offList); hold on;                                          %Plot to see the change in offdiagonal norm as iterations increase.
xlabel('Iterations'); ylabel('Off-diagnol Norm of A'); 

realEigenVal = eig(A);

eigenval = sort(eigenval)
realEigenVal = sort(realEigenVal)

toc;

