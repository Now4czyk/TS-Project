close all
clear

%% Parametry
m = .00575;
g = 9.81;
FemP1 = .017521;
FemP2 = .0058231;
f1 = .00014142;
f2 = .0045626;
ki = 2.5165;
ci = .0243;
d = .0792;
bd = .06;
xd = d-bd;
x10 = .008;
x20 = 0;
x30 = .75;
u10 = 0.2883767137;

a21 = 1/(2*m)*x30^2*FemP1/FemP2^2*exp(-x10/FemP2); 
a23 = -1/(m)*x30*FemP1/FemP2*exp(-x10/FemP2); 
a31 = 1/f1*exp(x10/f2)*(ki*u10+ci-x30); 
a33 = -f2/f1*exp(x10/f2);
b31 = ki*f2/f1*exp(x10/f2);
A = [0 1 0; a21 0 a23; a31 0 a33];
B = [0 0; 0 g; b31 0];
C = [1 0 0];
D = 0;

%% Dyskretyzacja modelu metodą Eulera w przód(na podstawie instrukcji lab5)
Tp = 0.001;
Ad = (eye(3)+A*Tp);
Bd = B*Tp;
Cd = C;
Dd = D;

