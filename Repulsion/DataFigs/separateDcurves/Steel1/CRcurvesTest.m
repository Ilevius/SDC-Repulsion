close all; clc;

simpleDC = load("Dcurves\simpleDcurves.txt");

dt1 = load("CDcurves\1.txt");
dt2 = load("CDcurves\2.txt");
dt3 = load("CDcurves\3.txt");
dt4 = load("CDcurves\4.txt");
dt5 = load("CDcurves\5.txt");
dt6 = load("CDcurves\6.txt");

IMMI2024(160 , 100, 14, 2, 7);
plot(dt1(:,3), dt1(:,1), '-b', dt1(:,3), dt1(:,2), '--b', ...
    dt2(:,3), dt2(:,1), '-m', dt2(:,3), dt2(:,2), '--m', dt3(:,3), dt3(:,1), '-g', dt3(:,3), dt3(:,2), '--g', ...
    dt4(:,3), dt4(:,1), '-r', dt4(:,3), dt4(:,2), '--r', dt5(:,3), dt5(:,1), '-c', dt5(:,3), dt5(:,2), '--c', ...
    dt6(:,3), dt6(:,1), '-', dt6(:,3), dt6(:,2), '--')
xlim([0 1.5]); ylim([-2 10]);
grid on;


% plot(simpleDC(:,1), simpleDC(:,2), '.', dt1(:,3), dt1(:,1), '-b', dt1(:,3), dt1(:,2), '--b', ...
%     dt2(:,3), dt2(:,1), '-y', dt2(:,3), dt2(:,2), '--y', dt3(:,3), dt3(:,1), '-g', dt3(:,3), dt3(:,2), '--g', ...
%     dt4(:,3), dt4(:,1), '-r', dt4(:,3), dt4(:,2), '--r', dt5(:,3), dt5(:,1), '-c', dt5(:,3), dt5(:,2), '--c', ...
%     dt6(:,3), dt6(:,1), '-', dt6(:,3), dt6(:,2), '--')


% plot(simpleDC(:,1), simpleDC(:,2), '.', dt4(:,3), dt4(:,1), '-r', dt4(:,3), dt4(:,2), '--r',...
%     dt5(:,3), dt5(:,1), '-g', dt5(:,3), dt5(:,2), '--g')
