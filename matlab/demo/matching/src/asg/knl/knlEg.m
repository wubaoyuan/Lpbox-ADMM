function [KP, KQ] = knlEg(XPs, XQs, egAlg)
% Compute edge affinity.
%
% Input
%   XPs     -  node feature, 1 x 2 (cell), dP x ni
%   XQs     -  edge feature, 1 x 2 (cell), dQ x mi
%   egAlg   -  method of computing edge kernel, 'toy' | 'cmum' | ...
%
% Output
%   KP      -  
%   KQ      -  edge affinity, m1 x m2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-10-2011

% dimension
[dEg, m1] = size(XEgs{1});
m2 = size(XEgs{2}, 2);

if strcmp(egAlg, 'toy')
    DEg = conDst(XEgs{1}, XEgs{2});
    KEg = exp(-DEg / .15);
    
elseif strcmp(egAlg, 'cmum')
    DEg = conDst(XEgs{1}, XEgs{2});
    KEg = exp(-DEg / 2500);
%    sigma = bandG(DEg, .2);    
%    KEg = exp(-DEg / 2 * sigma ^ 2);
    
elseif strcmp(egAlg, 'pas')
    
    % distance
    Dst1 = repmat(XEgs{1}(1, :)', 1, m2);
    Dst2 = repmat(XEgs{2}(1, :), m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(XEgs{1}(2, :)', 1, m2);
    Ang2 = repmat(XEgs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % default wEg
    if isempty(wEg)
        wEg = ones(dEg, 1);
    end
    
    % combine distance and angle
    KEg = exp(-wEg(1) * Dst - wEg(2) * Ang);
    
else
    error('unknown algorithm: %s', egAlg);
end
