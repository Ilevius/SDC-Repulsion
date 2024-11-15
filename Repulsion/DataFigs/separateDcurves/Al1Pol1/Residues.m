clc; close all; clear;

data1 = load("Residues\1.txt");
data2 = load("Residues\2.txt");
data3 = load("Residues\3.txt");
data4 = load("Residues\4.txt");
data5 = load("Residues\5.txt");
data6 = load("Residues\6.txt");
data7 = load("Residues\7.txt");
data8 = load("Residues\8.txt");
data9 = load("Residues\9.txt");
data10 = load("Residues\10.txt");



IMMIstyle2024(160 , 100, 14, 2, 7);
% IMMIstyle2024(widthmm , hightmm, textpt, lineWidth, markerSize)
marker = '-';
plot(data1(:,1), data1(:,3), marker,...
    data2(:,1), data2(:,3), marker,...
    data3(:,1), data3(:,3), marker, ...
    data4(:,1), data4(:,3), marker,...
    data5(:,1), data5(:,3), marker,...
    data6(:,1), data6(:,3), marker, ...
    data7(:,1), data7(:,3), marker,...
    data8(:,1), data8(:,3), marker,...
    data9(:,1), data9(:,3), marker, ...
    data10(:,1), data10(:,3), marker);
xlim([0 1.5]); ylim([0 0.05]);
grid on;





