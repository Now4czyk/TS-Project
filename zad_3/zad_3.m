close all; clear; clc;

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

%% definicja parametrów
a21 = 1/(2*m)*x30^2*FemP1/FemP2^2*exp(-x10/FemP2); 
a23 = -1/(m)*x30*FemP1/FemP2*exp(-x10/FemP2); 
a31 = 1/f1*exp(x10/f2)*(ki*u10+ci-x30); 
a33 = -f2/f1*exp(x10/f2);
b31 = ki*f2/f1*exp(x10/f2);

%% wyświetlenie transmitancji

A = [0 1 0; a21 0 a23; a31 0 a33];
B = [0 0; 0 g; b31 0];
C = [1 0 0];
D = 0;

ob = ss(A, B, C, D);
tf(ob)

%% Uruchomienie symulacji i stworzenie wykresów
out = sim("zad_3_sch.slx", .2);
pos = out.position.signals.values;
posTime = out.position.time;
reference = out.reference.signals.values;
referenceTime = out.reference.time;
control = out.control.signals.values;
controlTime = out.control.time;

plot(posTime, pos, referenceTime, reference);
xlabel("t [s]");
ylabel("position [m]");
legend('Sygnał wyjściowy', 'Sygnał referencyjny')
title("Sygnał wyjściowy (śledzenie sygnału sinusoidalnego)");
figure;
plot(controlTime, control);
xlabel("t [s]");
title("Sygnał sterujący (śledzenie sygnału sinusoidalnego)");
