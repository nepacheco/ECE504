clc; clear; close all;
%% Problem 1
A1 = [2 3 2; 3 1 0; 2 0 2]
eig(A1)

A2 = [0 0 -1; 0 0 0; -1 0 2]
eig(A2)

syms a1 a2 a3
A3 = [a1*a1 a1*a2 a1*a3; a2*a1 a2*a2 a2*a3; a3*a1 a3*a2 a3*a3]
simplify(eig(A3),'Steps', 5)

%% Problem 5

A = [-1 -1; 3 -5]
[Q, L] = eig(A)
% % syms a11 a12 a21 a22 real
% % syms m11 m12 m21 m22 real
% % syms u real
% % syms n11 n12 n21 n22 real
% % assume(m12 == m21)
% % assume(n12 == n21)
% % A = [a11 a12; a21 a22]
% % M = [m11 m12; m21 m22]
% % N = [n11 n12; n21 n22]
% % 
% % first = simplify(A'*M + M*A,'Steps', 10)
% % second = simplify(-N - 2*u*M,'Steps', 10)
% % eqn = A'*M + M*A + 2*u*M == -N;
% % solve(eqn, A)
