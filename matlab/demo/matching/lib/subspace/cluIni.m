function G0s = cluIni(K, para, varargin)
% Initalize the clustering.
%
% Input
%   K       -  frame similarity matrix, n x n
%   para    -  segmentation parameter
%     k     -  #clsuters
%     ini   -  initialization method, 'r' | 'sc'
%              'r': rand segmentation
%              'p': propagative segmentation
%     nIni  -  #initialization
%   varargin
%     save option
%     GT    -  ground-truth clustering, {[]}
%
% Output
%   G0s     -  initial clustering, 1 x nIni (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-12-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% save option
[savL, path] = psSave(varargin, 'subx', 'clu_ini');
if savL == 2 && exist(path, 'file')
    G0s = matFld(path, 'G0s');
    prom('m', 'old clu ini (%s): %d times\n', para.ini, para.nIni);
    return;
end
prom('m', 'new clu ini (%s): %d times\n', para.ini, para.nIni);

% function option
GT = ps(varargin, 'GT', []);

ini = para.ini; nIni = para.nIni; k = para.k;
G0s = cellss(1, nIni);

n = size(K, 1);
for i = 1 : nIni
    if strcmp(ini, 'r')
        G0s{i} = cluR(n, k);
    elseif strcmp(ini, 'sc')
        G0s{i} = cluSc(K, k);
    else
        error('unknown method');
    end

    % match with ground-truth segmentation
    if ~isempty(segT)
        G0s{i} = match(G0s{i}, 'GT', GT);
    end
end

% save
if savL > 0
    save(path, 'G0s');
end
