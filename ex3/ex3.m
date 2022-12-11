close all
% Trzeba wcześniej uruchomić symulację ex3PID.slx

pos = out.position.signals.values;
posTime = out.position.time;
reference = out.reference.signals.values;
referenceTime = out.reference.time;

plot(posTime, pos, referenceTime, reference)
legend("sygnał wyjściowy", "sygnał referencyjny")
xlabel("t [s]")
ylabel("position [cm]")