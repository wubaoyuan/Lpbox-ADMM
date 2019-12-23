function saveFigureAux(figureNumber,file);
% to make sure image does not resize while printing
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

set(figureNumber, 'PaperPositionMode', 'auto');

% otherwise white plots might turn to white
set(figureNumber,'InvertHardcopy','off');

figureColor = get(figureNumber,'Color');
set(figureNumber,'Color',[1,1,1]);

[pathName,fileName,ext] = fileparts(file);

try
    % dberror = disabledberror2;
    [dberror, dbstat] = disabledberror2;
catch
end

switch ext
    case '.eps'
        print( figureNumber, '-depsc', file ) ;
        %print( figureNumber, '-depsc','-tiff', file ) ;

        %         print(figureNumber, fullfile(pathName,[fileName,'.jpg']) , '-djpeg', '-r0');
        %         print(figureNumber, fullfile(pathName,[fileName,'.png']) , '-dpng', '-r0');
        print(figureNumber, fullfile(pathName,[fileName,'.png']) , '-dpng', '-r300');
    case ''
        ext = '.jpg';
        saveas(figureNumber,fileName);
    case '.png'
        print(figureNumber, file , '-dpng', '-r0');
    case {'.jpg','.jpeg'}
        print(figureNumber, file , '-djpeg', '-r0');
    otherwise
        %     print(figureNumber, file , '-dpng', '-r0'  ,  '-loose');
        %             print(figureNumber, file , '-dpng', '-r0');
        %     print(figureNumber, file , '-depsc', '-r0',  '-loose');
        saveas(figureNumber,file);
end
try
    enabledberror2(dberror, dbstat);
catch
end
set(figureNumber,'Color',figureColor);
