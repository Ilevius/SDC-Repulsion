function IMMIstyle2024(widthmm , hightmm, textpt, lineWidth, markerSize)
multiple = 1;
widpi = (650/127)*widthmm*multiple;
higpi = (480/127)*hightmm*multiple;
textri = textpt*multiple;
set(0, 'DefaultAxesFontSize', textri, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', textri, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', lineWidth);

set(gcf,'color','w');
set(gca, 'GridAlpha', 0.2);
set(groot,'defaultLineMarkerSize',markerSize);


aPlot = figure();
set (aPlot, 'Position', [200 200 widpi higpi]);

end