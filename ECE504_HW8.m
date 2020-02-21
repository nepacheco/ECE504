clc; clear;close all

%% Problem 1
A = [-2 3 0; 1 0 0; 0 1 0];
B = [1;0;0];
C = [0 1 3];

control = ctrb(A,B)
disp("rank(control)")
rank(control)

syms lambda real
control2 = [A-lambda*eye(3,3) B];
disp("rank(control)")
rank(control2)

observation = obsv(A,C)
disp("rank(observation)")
rank(observation)

%% Problem 2

A = [0 1 0; 0 0 1; 0 2 -1];
B = [0 2;1 0;0 0];
C = [1 0 -2];

control = ctrb(A,B)
rank(control)

syms lambda real
control2 = [A-lambda*eye(3,3) B]
rank(control2)

observation = obsv(A, C);
rank(observation)

%% Problem 4

syms A11 A12 A21 A22 B1
A = [A11 A12;A21 A22];
B = [B1;0];

control = [B A*B]
other  =[A22 A22*A21]