function cnames = className(cnames0, CMap)
% Pick out the used class names, and fill them in a new name list. 
%
% Input
%   cnames0  -  original class names, 1 x maxc (cell)
%   CMap     -  dictionary, k x maxc
%
% Output
%   cnames   -  new class names, 1 x k (cell)
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

k = size(CMap, 1);
cnames = cell(1, k);
for c = 1 : k
    cnames{c} = cnames0{find(CMap(c, :))};
end
