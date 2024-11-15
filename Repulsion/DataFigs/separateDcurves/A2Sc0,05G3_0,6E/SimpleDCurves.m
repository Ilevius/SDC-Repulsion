clc; close all; clear;

simpleDC = load("Dcurves\simpleDcurves.txt");
plot(simpleDC(:,1), simpleDC(:,2), '.')
