close all; clc;

simpleDC = load("Dcurves\simpleDcurves.txt");

dt = load("Dcurves\CRcurve.txt");


plot(simpleDC(:,1), simpleDC(:,2), '.')
hold on
plot(dt(:,3), dt(:,1), '.', dt(:,3), dt(:,2), 'x')