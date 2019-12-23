function CMUM = cmumHuman
% Load ground label file for CMU Motion sequence.
%
% Output
%   CMUM   -  a container
%     tag  -  sequence name
%       seg
%       cnames
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-03-2013

% specified in addPath.m
global footpath;
foldpath = sprintf('%s/data/cmum', footpath);
matpath = sprintf('%s/cmum.mat', foldpath);

% if mat existed, just load
if exist(matpath, 'file')
    CMUM = matFld(matpath, 'CMUM');
    prInOut('cmumHuman', 'old');
    return;
end
prIn('cmumHuman', 'new');

% all subject
tags = getAllNames;

m = length(tags);

for i = 1 : m
    tag = tags{i};
    
    if strcmp(tag, 'house')
        matpathi = sprintf('%s/%s/label/st.mat', foldpath, tag);
        Tmp = load(matpathi);
        P = Tmp.P;
    else
        P = cmumLabel(tag);
    end
    
    
%    equal('P', Tmp.P, P0s);
    
    % markers' position
    Xis = P;
    
    for i = 1 : length(Xis)
        Xis{i} = Xis{i}';
    end
    
    % store
    CMUM.(tag).XTs = Xis;
end

% save
save(matpath, 'CMUM');

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tags = getAllNames
% Obtain all names under the folder data/cmum.
%
% Output
%   tags  -  names, 1 x m (cell)

tags = {'house', 'hotel'};
%tags = {'house'};
