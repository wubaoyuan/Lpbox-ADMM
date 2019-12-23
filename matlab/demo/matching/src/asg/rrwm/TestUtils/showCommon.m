if ~exist('nFigureHandle')
    nFigureHandle = 0;
end
nFigureHandle = nFigureHandle + 1;
figureHandle(nFigureHandle) = figure('Position', figureSize);

hold on
for k = 1:length(settings.Algorithm)
    plot(settings.Con.var, meanData(:,k), 'color', settings.Algorithm(k).color, ...
        'linestyle', settings.Algorithm(k).lineStyle, 'linewidth', 3, ...
        'Marker', settings.Algorithm(k).Marker, 'MarkerSize', sizeMarker);
end

xlabel(settings.Con.str, 'fontsize', sizeFont);
ylabel(['\fontname{times new roman}' yLabelText], 'fontsize', sizeFont);

Xmin = min(settings.Con.var); Xmax = max(settings.Con.var);
% Ymin = min(meanData(:)); Ymax = max(meanData(:));
axis([Xmin Xmax Ymin Ymax*1.01]);

if bDisplayLegend
    h = legend(settings.Algorithm(:).strName, 'location', 'best');
    set(h, 'fontsize', legendFontSize);
end

for k = length(settings.Algorithm):-1:1
    plot(settings.Con.var, meanData(:,k), 'color', settings.Algorithm(k).color, ...
        'linestyle', settings.Algorithm(k).lineStyle, 'linewidth', 3, ...
        'Marker', settings.Algorithm(k).Marker, 'MarkerSize', sizeMarker);
end

for k = 1:length(settings.Fix)
    text(0.1*(Xmax-Xmin)+Xmin, stepSize*(k)*(Ymax-Ymin)+Ymin, ...
        ['\fontname{times new roman}' settings.Fix(k).str '\fontname{helvetica}\rm = ' num2str(settings.Fix(k).var)], 'fontsize', sizeFont);
end
text(0.1*(Xmax-Xmin)+Xmin, stepSize*(k+1)*(Ymax-Ymin)+Ymin, ...
    ['\fontname{times new roman}# of inlier \it{n_{in}}\fontname{helvetica}\rm = ' num2str(settings.nInlier)], 'fontsize', sizeFont);

hold off