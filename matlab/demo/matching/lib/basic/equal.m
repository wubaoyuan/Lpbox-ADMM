function flag = equal(s, a, b, varargin)
% Testify whether two objs are equal.
% The obj can be 'numeric', 'cell', 'char' and 'struct'.
%
% Input
%   s       -  name
%   a       -  1st obj
%   b       -  2nd obj
%   varargin
%     pr    -  prompt flag, {'y'} | 'n'
%
% Output
%   flag    -  boolean flag
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-13-2012

% function option
isPr = psY(varargin, 'pr', 'y');

[flag, dif] = eqAB(a, b);
if isPr
    if flag
        fprintf('%s equal\n', s);
    else
        fprintf('%s inequal, %s\n', s, dif);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag, dif] = eqAB(a, b)

% find the type
types = {'numeric', 'cell', 'char', 'struct'};
k = length(types);

for c = 1 : k
    isA = isa(a, types{c});
    isB = isa(b, types{c});
    
    if isA && isB
        break;
    end

    if isA ~= isB
        flag = false;
        dif = 'different type';
        return;
    end
end

% numeric matrix
if strcmp(types{c}, 'numeric')
    % dimension
    flag = eqDim(a, b);
    if ~flag
        dif = 'different dimension';
        return;
    end
    
    % infinite value
    a = a(:);
    b = b(:);
    visa = isinf(a);
    visb = isinf(b);
    flag = all(visa == visb);
    if ~flag
        dif = 'different (in)finite';
        return;
    end
    a = a(~visa);
    b = b(~visb);
    
    % nan value
    a = a(:);
    b = b(:);
    visa = isnan(a);
    visb = isnan(b);
    flag = all(visa == visb);
    if ~flag
        dif = 'different (in)nan';
        return;
    end
    a = a(~visa);
    b = b(~visb);

    if isempty(a)
        flag = true;
    else
        d = max(abs(a - b));
        flag = d < 1e-7;
        if ~flag
            dif = sprintf('different value %.2e', d);
            return;
        end
    end

% cell
elseif strcmp(types{c}, 'cell')
    % dimension
    flag = eqDim(a, b);
    if ~flag
        dif = 'different dimension';
        return;
    end

    % content
    flag = true;
    for i = 1 : numel(a)
        if ~eqAB(a{i}, b{i})
            flag = false;
            dif = 'different value';
            return;
        end
    end

% char
elseif strcmp(types{c}, 'char')
    flag = strcmp(a, b);
    if ~flag
        dif = sprintf('different string: %s vs %s', a, b);
        return;
    end

% struct array
elseif strcmp(types{c}, 'struct')
    % dimension
    flag = eqDim(a, b);
    if ~flag
        dif = 'different dimension';
        return;
    end
    
    % fieldnames
    nmAs = fieldnames(a);
    nmBs = fieldnames(b);
    flag = length(nmAs) == length(nmBs);
    if ~flag
        dif = 'different number of fields';
        return;
    end

    % content
    for j = 1 : length(nmAs)
        nmA = nmAs{j};
        nmB = nmBs{j};
        flag = strcmp(nmA, nmB);
        if ~flag
            dif = sprintf('different field name: %s vs %s', nmA, nmB);
            return;
        end
    
        for i = 1 : numel(a)
            if ~eqAB(a(i).(nmA), b(i).(nmB))
                flag = false;
                dif = sprintf('different filed value for field: %s', nmA);
                return;
            end
        end
    end
end

dif = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag = eqDim(a, b)

siza = size(a);
sizb = size(b);

if length(siza) ~= length(sizb)
    flag = false;
else
    flag = all(siza == sizb);
end
