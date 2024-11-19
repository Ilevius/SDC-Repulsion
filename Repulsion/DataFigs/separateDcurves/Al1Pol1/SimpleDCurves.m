clc; close all; clear;

RPDots = load("Dcurves/simpleDcurves.txt");
RPole1 = load("Dcurves/1.txt");
RPole2 = load("Dcurves/2.txt");
RPole3 = load("Dcurves/3.txt");
plot(RPDots(:,1), RPDots(:,2), '.', RPole1(:,1), RPole1(:,2), RPole2(:,1), RPole2(:,2), RPole3(:,1), RPole3(:,2))

