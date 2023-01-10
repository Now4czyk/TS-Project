close all; clear;

%% Parametry
m = .00575;
g = 9.81;
FemP1 = .017521;
FemP2 = .0058231;
f1 = .00014142;
f2 = .0045626;
ki = 2.5165;
ci = .0243 * sign(4-2);
d = .0792;
bd = .06;
xd = d-bd;
x10 = .008;
x20 = 0;
x30 = .75;
u10 = 1/ki*(x30-ci);

a21 = 1/(2*m)*x30^2*FemP1/FemP2^2*exp(-x10/FemP2); 
a23 = -1/(m)*x30*FemP1/FemP2*exp(-x10/FemP2); 
a31 = 1/f1*exp(x10/f2)*(ki*u10+ci-x30); 
a33 = -f2/f1*exp(x10/f2);
b31 = ki*f2/f1*exp(x10/f2);
A = [0 1 0; a21 0 a23; a31 0 a33];
B = [0 0; 0 g; b31 0];
C = [1 0 0];
D = 0;

%% Postać kanoniczna sterowalna
polyA = -poly(A);
A2 = [0 1 0; 0 0 1; polyA(4) polyA(3) polyA(2)];
B2 = [0; 0; 1];

% Macierz podobieństwa
S1 = ctrb(A,B(:, 1));
S2 = ctrb(A2, B2);
P=S2*S1^(-1);

alpha = P*A*P^(-1)
beta = P*B(:, 1);
beta(:, 2) = [0; 0; g];
beta
gamma = C*P^(-1)
delta = D

%% Postać diagonalna
[M, b] = jordan(A);

P = M^(-1);

alpha = P*A*P^(-1)
beta = P*B
gamma = C*P^(-1)
delta = 0

