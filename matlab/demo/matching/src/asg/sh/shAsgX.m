function shAsgX(asgs, Ax, algs)
% Show the correspondence matrix.
% 
% Input
%   asgs    -  assignment, 1 x mAlg (cell)
%   Ax      -  axis, rows x cols (cell)
%   algs    -  algorithm name
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-26-2012

% dimension
mAlg = length(asgs);

for iAlg = 1 : mAlg
    asg = asgs{iAlg};
    shM(asg.X, 'ax', Ax{iAlg});
    
    if ~exist('algs', 'var') || isempty(algs)
        alg = asg.alg;
    else
        alg = algs{iAlg};
    end
    
    if isfield(asg, 'acc')
        if isfield(asg, 'obj')
            title(sprintf('%s: acc %.2f obj %.2f', alg, asg.acc, asg.obj));
        else
            title(sprintf('%as: acc %.2f', alg, asg.acc));            
        end
    else
        title(sprintf('%s', alg));
    end
end

% close other axes
[rows, cols] = size(Ax);
for i = mAlg + 1 : rows * cols
    set(Ax{i}, 'Visible', 'off');
end