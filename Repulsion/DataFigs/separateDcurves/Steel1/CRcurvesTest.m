close all; clc;

simpleDC = load("Dcurves\simpleDcurves.txt");

dt1 = load("CDcurves\1.txt");
dt2 = load("CDcurves\2.txt");
dt3 = load("CDcurves\3.txt");
dt4 = load("CDcurves\4.txt");


plot(simpleDC(:,1), simpleDC(:,2), '.', dt1(:,3), dt1(:,1), '-', dt1(:,3), dt1(:,2), '--', ...
    dt2(:,3), dt2(:,1), '-', dt2(:,3), dt2(:,2), '--', dt3(:,3), dt3(:,1), '-', dt3(:,3), dt3(:,2), '--', ...
    dt4(:,3), dt4(:,1), '-', dt4(:,3), dt4(:,2), '--')
