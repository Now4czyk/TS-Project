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

%% Sprawdzenie stabilności, sterowalności i obserwowalności
a21 = 1/(2*m)*x30^2*FemP1/FemP2^2*exp(-x10/FemP2); 
a23 = -1/(m)*x30*FemP1/FemP2*exp(-x10/FemP2); 
a31 = 1/f1*exp(x10/f2)*(ki*u10+ci-x30); 
a33 = -f2/f1*exp(x10/f2);
b31 = ki*f2/f1*exp(x10/f2);

%% Uruchomienie symulacji i stworzenie wykresów

out = sim("ex3PID.slx", 40);

pos = out.position.signals.values;
posTime = out.position.time;
reference = out.reference.signals.values;
referenceTime = out.reference.time;

plot(posTime, pos, referenceTime, reference)
legend("sygnał wyjściowy", "sygnał referencyjny")
xlabel("t [s]")
ylabel("position [cm]")