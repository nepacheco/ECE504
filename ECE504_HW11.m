clc; clear; close all

%% Problem 1
A = [ 1 1 -2; 0 1 1; 0 0 1];
B = [1; 0; 1];
C = [2 0 0 ];

syms L
act_poly = simplify(det(L*eye(3,3) - A), 'Steps', 10,'All', true)
Abar = [3 -3 1; 1 0 0; 0 1 0];
Bbar = [1; 0 ; 0];

control_mat = ctrb(A,B)
ccf_cont_mat = ctrb(Abar, Bbar)
P = ccf_cont_mat * pinv(control_mat)

des_poly = conv(conv([1 2],[1 (1 -i)]),[1 1+i])

G = [4+3 6-3 4+1]
K = G*P
% Verification
p = [-2 -1+i -1-i]
place(A,B,p)

syms p
As = A-B*K;
eigs = eig(As)
tf = C*inv(L*eye(3,3) - As)*p*B

%% Problem 2
A = [1 1 -2 ; 0 1 1 ; 0 0 1];
B = [1; 0 ; 1];
C = [2 0 0];

syms L
act_poly = simplify(det(L*eye(3,3) - A), 'Steps', 10,'All', true)
Abar = [3 -3 1; 1 0 0; 0 1 0]
Bbar = [1; 0 ; 0]

control_mat = ctrb(A,B)
ccf_cont_mat = ctrb(Abar, Bbar)
P = ccf_cont_mat * inv(control_mat)

des_poly = conv(conv([1 0],[1 0]),[1 0])

G = [3 -3 1]
K = G*P
As = A-B*K
eig(As);

syms x1 x2 x3 real
x_0 = [x1; x2; x3];
x_1 = As*x_0
x_2 = As*x_1
x_3 = As*x_2


%% Problem 3

A = [1 1 -2 ; 0 1 1 ; 0 0 1];
B = [1; 0 ; 1];
C = [2 0 0];
K = [1 5 2];
syms z
As = (A-B*K);
tf = C*inv(z*eye(3,3)-As)*B;
simplify(tf,'Steps',10,'All',true)

syms r x1 x2 x3 real
x = [];
x = [x1;x2;x3];
y =[];
p = 0.5;
for k = [0,1,2,3,4]
    disp('K value')
    disp(k);
    n = k + 1;
    x = [x As*x(:,n) + p*B*r]
    y = [y C*x(:,n)]
end

%% Problem 4
A = [2 1 0 0 ; 0 2 0 0 ; 0 0 -1 0 ; 0 0 0 -1];
B = [0 ; 1 ; 1 ; 1];
syms L
act_poly = simplify(det(L*eye(4,4) - A), 'Steps', 10,'All', true)

Accf = [ 2 3 -4 -4; 1 0 0 0; 0 1 0 0 ; 0 0 1 0 ];
Bccf = [1 ; 0 ; 0 ; 0];
ccf_contr_mat = ctrb(Abar,Bbar)
contr_mat = ctrb(A,B)
li_columns = rank(contr_mat);
Q1 = contr_mat(:,1:li_columns);
Q2 = null(contr_mat(:,1:li_columns)');

Q = [Q1 Q2];
P = inv(Q);
Abar = P*A*Q;
Bbar = P*B;

Abarc = Abar(1:3,1:3);
A12 = Abar(1,4);
Abar_nc = Abar(4,4);
Bbarc = Bbar(1:3,1);
Bbar_nc = Bbar(4,1);

eig(Abar_nc)
eig(Abarc)


%% Problem 5
A = [2 1 ; -1 1];
B = [1 ; 2];
C = [1 1];
syms L 
tf = simplify(det(L*eye(2,2) - A'),'Steps', 10)
contr_mat = ctrb(A',C')
rank(contr_mat)

Abar = [3 -3; 1 0]
Bbar = [1 ; 0]

ccf_cont_mat = ctrb(Abar, Bbar)

P = ccf_cont_mat * inv(contr_mat);
desired_tf = conv([1 2-2*i], [1 2+2*i])
G = [4+3 8-3];
K = G*P;
L = K'

As = A-L*C
eig(As)

%% Problem 5 Lyaponov Method
A = [2 1 ; -1 1];
B = [1 ; 2];
C = [1 1];
eig(A)
eig([-2 2; -2 -2])

mat = [4 -1 -2 0 1;
       1 3 0 -2 1;
       2 0 4 -1 0;
       0 2 1 3 0]
soln = rref(mat);
T = [soln(1,5) soln(2,5); soln(3,5) soln(4,5)];
F = [-2 2; -2 -2];
l = [1;0];

inv(T)*F*T
inv(T)*T*B
inv(T)*l
