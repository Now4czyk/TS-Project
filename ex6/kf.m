function x = kf(we)
%% Odczyt wejsc
global P
global P_post
global x_post
y = we(1);  % aktualny pomiar
u = we(2);  % aktualne wymuszenie
x1 = we(3:5);   % poprzedni wektor stanu
yy = we(6);

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

%% Dyskretyzacja modelu metodą Eulera w przód(na podstawie instrukcji lab5)
Tp = 0.001;
Ad = (eye(3)+A*Tp);
Bd = B*Tp;
Cd = C;
Dd = D;


P_post = I*1000;
x_post = [xd; 0; 0];

%% Nastawy filtru
var_n = % wariancja szumu pomiarowego
var_v = % wariancje szumow wewnetrznych

Q = diag(var_v);    % macierz kowariancji szumow wewnetrznych
R = diag(var_n);  	% macierz kowariancji szumow pomiarowych


%% Algorytm filtru Kalmana
% Etap predykcji
x_d = Ad*x_post + Bd*u;
P_d = Ad*P_post*Ad' + Q;

% Etap filtracji
K = P_d*Cd' * (Cd*P_d*Cd' + R)^-1;
x_hat = x_d + K*(yy - Cd*x_d);
P = (I - K*Cd)*P_d;
    
% Przekazanie danych
x_post = x_hat
P_post = P;

x = x_hat
end

