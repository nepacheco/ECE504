%% ECE Homework 9
clc; clear; close all

%% Problem 1
A = [-1 4; 4 -1];
B = [1;1];
C = [1 1];

[ABAR,BBAR,CBAR,T,K] = ctrbf(A,B,C);

%% Problem 2

syms L1 L2 real

A = [L1 1 0 0 0;
     0 L1 1 0 0;
     0 0 L1 0 0;
     0 0 0 L2 1;
     0 0 0 0 L2];
B = [0;1;0;0;1];
C = [0 1 1 0 1];

% initial controllability matrix
cont = simplify([B A*B A*A*B A*A*A*B A*A*A*A*B], 'Steps', 10);
disp("Rank of original Controllability Matrix")
rank(cont)

% Creating Q1 and Q2 from controllability matrix
Q1 = [0,  1, 2*L1, 3*L1^2;
    1, L1, L1^2,   L1^3;
    0,  0,    0,      0;
    0,  1, 2*L2, 3*L2^2;
    1, L2, L2^2,   L2^3];
Q2 = [0; 0; 1; 0; 0];

% disp("Q matrix");
Q = [Q1 Q2];
disp("Rank Q");
rank(Q)

% disp("P matrix");
P = inv(Q);

Abar = simplify(P*A*Q, 'Steps', 10);
Bbar = simplify(P*B, 'Steps', 10);
Cbar = simplify(C*Q, 'Steps', 10);

Abar_c = [ 0, 0, 0,              -L1^2*L2^2;
           1, 0, 0,       2*L1*L2*(L1 + L2);
           0, 1, 0, - L1^2 - 4*L1*L2 - L2^2;
           0, 0, 1,             2*L1 + 2*L2];
       
Bbar_c = [1;0;0;0];
Cbar_c = [ 2, L1 + L2, L1^2 + L2^2, L1^3 + L2^3];

disp("Rank of controllability matrix after controllability decomposition")
control_decomp = [Bbar_c Abar_c*Bbar_c Abar_c*Abar_c*Bbar_c...
                    Abar_c*Abar_c*Abar_c*Bbar_c];
rank(control_decomp)

% observability matrix after doing controllability decomposition
observ = simplify([Cbar_c;
    Cbar_c*Abar_c;
    Cbar_c*Abar_c*Abar_c;
    Cbar_c*Abar_c*Abar_c*Abar_c], 'Steps', 10);
disp("rank of observability matrix after doing controllability decomposition")
rank(observ)

P1_c = [           2,     L1 + L2, L1^2 + L2^2, L1^3 + L2^3
          L1 + L2, L1^2 + L2^2, L1^3 + L2^3, L1^4 + L2^4];

P2_c = [L1*L2, -L1-L2, 1 0;
        L1^2*L2 + L1*L2^2, -L1^2 - L1*L2 - L2^2, 0 1];
P_c = [P1_c; P2_c];
disp("rank of P")
rank(P_c)
% disp("Comparing Q_c with P1_c inverse")

Q_c = simplify(inv(P_c), 'Steps', 100);
Q1_c = simplify(pinv(P1_c), 'Steps', 100);

Abar2 = simplify(P_c*Abar_c*Q_c, 'Steps', 10);
Bbar2 = simplify(P_c*Bbar_c, 'Steps', 100);
Cbar2 = simplify(Cbar_c*Q_c, 'Steps', 100);

Abar2_o = simplify(P1_c*Abar_c*Q_c(:,1:2), 'Steps', 10)
Bbar2_o = simplify(P1_c*Bbar_c, 'Steps', 100)
Cbar2_o = simplify(Cbar_c*Q_c(:,1:2), 'Steps', 100)

syms s
disp("Original Transfer Function")
orig_trans = simplify(C*inv(s*eye(5,5) - A)*B, 'Steps', 10)
disp("Kalman decomposition Transfer Function")
new_trans = simplify(Cbar2_o*inv(s*eye(2,2) - Abar2_o)*Bbar2_o, 'Steps', 10)
