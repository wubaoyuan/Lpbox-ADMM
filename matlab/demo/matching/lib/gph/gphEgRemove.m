function Eg = gphEgRemove(Eg0, X)
% Connecting points to obtain edge.
%
% Input
%   EgTr    -  edge of training graph, 2 x nTr
%   Ct      -  constraint, nTr x nTe
%
% Output
%   EgTe    -  edge of testing graph, 2 x nTe
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-12-2012

% dimension
idx = find(asg.X(:));
[~, iTes] = ind2sub(size(asg.X), idx);
A = repmat(Eg(1, :)', 1, 29);
B = repmat(iTes, length(Eg(1, :)), 1);

[nTr, nTe] = size(Ct);
mTr = size(EgTr, 2);
mS = round(nTe / nTr);

% index
IdxTe = zeros(nTr, mS);
for i = 1 : nTr
    IdxTe(i, :) = find(Ct(i, :) == 1);
end

Is = zeros(mS, mS, mTr);
Js = zeros(mS, mS, mTr);
for i = 1 : mTr
    x = EgTr(1, i);
    y = EgTr(2, i);
    
    idxIs = IdxTe(x, :);
    idxJs = IdxTe(y, :);
    
    [Is(:, :, i), Js(:, :, i)] = meshgrid(idxIs, idxJs);
end

EgTe = [Is(:), Js(:)]';
