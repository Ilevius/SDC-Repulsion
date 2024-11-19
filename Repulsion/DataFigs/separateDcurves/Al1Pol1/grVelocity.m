clc; close all; clear;

data1 = load("Dcurves\1.txt");
data2 = load("Dcurves\2.txt");
data3 = load("Dcurves\3.txt");
data4 = load("Dcurves\4.txt");
data5 = load("Dcurves\5.txt");
data6 = load("Dcurves\6.txt");

h = 0.01;  %h=0.005
[gv12, gv11] = dxdf(data1(:,1), data1(:,2), h);
[gv22, gv21] = dxdf(data2(:,1), data2(:,2), h);
[gv32, gv31] = dxdf(data3(:,1), data3(:,2), h);
[gv42, gv41] = dxdf(data4(:,1), data4(:,2), h);
[gv52, gv51] = dxdf(data5(:,1), data5(:,2), h);
[gv62, gv61] = dxdf(data6(:,1), data6(:,2), h);



mrk = '-'; yFac = 2*pi;
plot(gv11, gv12*yFac, mrk, gv21, gv22*yFac, mrk, gv31, gv32*yFac, mrk, gv41, gv42*yFac, mrk, ...
    gv51, gv52*yFac, mrk, gv61, gv62*yFac, mrk);