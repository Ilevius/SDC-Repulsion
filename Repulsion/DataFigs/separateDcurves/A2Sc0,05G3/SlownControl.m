clc; close all; clear;

data1 = load("Residues\1.txt");
data2 = load("Residues\2.txt");
data3 = load("Residues\3.txt");
data4 = load("Residues\4.txt");
data5 = load("Residues\5.txt");
data6 = load("Residues\6.txt");
data7 = load("Residues\7.txt");

simpleDC = load("Dcurves\simpleDcurves.txt");


marker = '-';
IMMIstyle2024(160 , 100, 14, 2, 7);
plot(data1(:,1), data1(:,2)./(data1(:,1)*2*pi), marker, data2(:,1), data2(:,2)./(data2(:,1)*2*pi), marker, data3(:,1), data3(:,2)./(data3(:,1)*2*pi), marker, ...
    data4(:,1), data4(:,2)./(data4(:,1)*2*pi), marker, data5(:,1), data5(:,2)./(data5(:,1)*2*pi), marker, data6(:,1), data6(:,2)./(data6(:,1)*2*pi), marker, ...
    data7(:,1), data7(:,2)./(data7(:,1)*2*pi), marker);
xlim([0 1.5]); ylim([0 0.5]);
grid on;

hold on;
plot(simpleDC(:,1), simpleDC(:,2)./simpleDC(:,1)/2/pi, '.');
