%% Set Varying Variables
sizeFont = 20;
bDisplayLegend = 1;
legendFontSize = 14;
sizeMarker = 8;
stepSize = 0.1;
figureSize = [1920*0.1 1200*0.1 1920*0.4 1200*0.4];

%%
meanData = meanScore;
Ymin = 0; Ymax = 1;
yLabelText = 'Accuracy';
showCommon;
%%
meanData = meanEnergy;
Ymin = min(meanData(:)); Ymax = max(meanData(:));
yLabelText = 'Objective score';
showCommon;
%%
meanData = meanTime;
Ymin = 0; Ymax = max(meanData(:));
yLabelText = 'Time';
showCommon;

clear Xmax Xmin Ymax Ymin figureSize legendFontSize nFigureHandle sizeFont sizeMarker stepSize yLabelText h meanData bDisplayLegend