clc; close all; clear;

data1 = load("Dcurves\1.txt");
data2 = load("Dcurves\2.txt");
data3 = load("Dcurves\3.txt");
data4 = load("Dcurves\4.txt");
data5 = load("Dcurves\5.txt");
data6 = load("Dcurves\6.txt");
data7 = load("Dcurves\7.txt");
data8 = load("Dcurves\8.txt");
data9 = load("Dcurves\9.txt");
data10 = load("Dcurves\10.txt");


marker = '-';
IMMIstyle2024(160 , 100, 14, 2, 7);
plot(data1(:,1), data1(:,2)./(data1(:,1)*2*pi), marker,...
    data2(:,1), data2(:,2)./(data2(:,1)*2*pi), marker,...
    data3(:,1), data3(:,2)./(data3(:,1)*2*pi), marker, ...
    data4(:,1), data4(:,2)./(data4(:,1)*2*pi), marker,...
    data5(:,1), data5(:,2)./(data5(:,1)*2*pi), marker,...
    data6(:,1), data6(:,2)./(data6(:,1)*2*pi), marker, ...
    data7(:,1), data7(:,2)./(data7(:,1)*2*pi), marker,...
    data8(:,1), data8(:,2)./(data8(:,1)*2*pi), marker,...
    data9(:,1), data9(:,2)./(data9(:,1)*2*pi), marker,...
    data10(:,1), data10(:,2)./(data10(:,1)*2*pi), marker);
% xlim([0 1.6]); ylim([0 0.5]);
grid on;


