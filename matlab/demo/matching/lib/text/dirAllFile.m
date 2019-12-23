function nms = dirAllFile(dpath, dpatt, fpatt, prex)
% Recursively parsing a directory to obtain the file tree.
%
% Input
%   dpath   -  path of a directory
%   dpatt   -  pattern of matched directory names
%   fpatt   -  pattern of matched file names
%   prex    -  prefix name, [] | string array
%
% Output
%   nms     -  list of file names, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% sub-folder
ds = dir(dpath);
n = length(ds);

% per sub-folder
nmss = cell(1, n);
vis = zeros(1, n);
for i = 1 : n
    dnm = lower(ds(i).name);
    
    % prefix
    if isempty(prex)
        prexi = dnm;
    else
        prexi = [prex '_' dnm];
    end
    
    % directory
    if ds(i).isdir
        pos = regexp(dnm, dpatt, 'once');
        if isempty(pos)
            continue;
        end
        
        nmss{i} = dirAllFile([dpath '/' dnm], dpatt, fpatt, prexi);
        
    % file
    else
        pos = regexp(dnm, fpatt, 'once');
        if isempty(pos)
            continue;
        end
        
        % remove .xxx
        pos = find(prexi == '.', 1, 'first');
        if ~isempty(pos)
            prexi(pos : end) = [];
        end

        nmss{i} = prexi;
    end
    vis(i) = 1;
end

% concatenate
nms = cellCat(nmss(vis == 1));
