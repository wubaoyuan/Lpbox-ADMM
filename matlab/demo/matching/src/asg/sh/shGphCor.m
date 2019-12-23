function ha = shGphCor(gphs, X, XT, parCor)
% Show correspondence.
%
% Input
%   gphs    -  graph, 1 x 2 (cell)
%   X       -  correspondence, n1 x n2
%   XT      -  ground-truth correspondence, [] | n1 x n2
%   parCor  -  parameter
%     cor   -  method for visualizing correspondence, {'ln'} | 'col'
%                'ln': using line
%                'col': using color
%
% Output
%   ha      -  handle
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-30-2012

% function option
cor = ps(parCor, 'cor', 'ln');

if strcmp(cor, 'ln')
    ha = shGphCorLn(gphs, X, XT, parCor);
elseif strcmp(cor, 'lnC')
    ha = shGphCorLnC(gphs, X, parCor);    
elseif strcmp(cor, 'col')
    ha = shGphCorCol(gphs, X, XT, parCor);
else
    error('unknown method: %s', cor);
end 
