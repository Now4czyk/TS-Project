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

%% Sprawdzenie jakie wartości wymuszeń spowodują przyciągnięcie kuli. Zbadanie odpowiedzi skokowej w układzie otwartym.
ob = ss(A, B, C, D);
sys1 = tf(ob)
% figure; step(sys1)

%% Wzmocnienia sprzężenia od stanu
% dla s = [2 7 80] odpowiedź ustala się na wartości
s = [600, 620, 630]';
k = -place(A, B(:,1), -s)
% A2=A+B*k
% eig(A2)
% figure; step(ss(A2, B, C, D));

%% liczenie stałowartościowej
% 
% u0 = .05;
% 
% % to na kartce, bo od razu wychodzi, nawet nie trzeba liczyć, pogczamp
% % super aleluja
x1 = .01;
x2 = 0;
x3 = -(a21*x1+g)/a23
% 
u0 = -(a31*x1+a33*x3)/b31;
% 
% % liczyłem to po północy, mogą być błędy
% % x3 = (g*a31-b31*a21)/(a33*a21-a23*a23)*u0
% % x1 = -(a23*x3)/a21
% % x2=0
% 
% xr = [x1, x2, x3]

%% Sprzężenie od stanu

%pochodne porównać do zera
%na u0+k wzorekr z kartkowki  kx-x=x0

% x0 = [0, 0, 0]   % warunki poczatkowe
% t = 0:0.001:30;     % wektor czasu
% [t, x] = ode45(@odefun, t, x0, [], A, B, k, u0, xr);
% 
% % Rysowanie wykresow
% figure; plot(t, x(:,1));
% xlabel('t [sek]'); ylabel('y(t) [cm]');
% 
% function dxdt = odefun(t, x, A, B, k, u0, xr)
%     u = u0 + k*(x-xr)
%     dxdt = A * x + B(:,1) * u;
% end

