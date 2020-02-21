syms a1 a2 a3 a4 c4 real
syms s
A = [-2 1 2 a1; 1 0 0 a2; 0 1 0 a3; 0 0 0 a4];
B = [1 ; 0 ; 0 ; 0];
C = [ 0 1 2 c4];

obs = [C; C*A; C*A*A; C*A*A*A];
rank(obs)

tf = C*inv(s*eye(4,4)- A)*B;
tf = simplify(tf, 'Steps', 10)