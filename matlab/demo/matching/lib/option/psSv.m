function [svL, path, prex, fold, subx] = psSv(option, varargin)
% Parse the save option and generate the output path.
% The input have 5 parameters:
%    'svL': save level, 0 | 1 | 2
%           0: do the algorithm every time and don't save
%           1: do the algorithm every time but save the result
%           2: do the algorithm only if the save file doesn't exist and if redo, save the result
%    'prex': prefix
%    'fold': fold name
%    'subx': subfix
%    'path': path (if exist, ignore the prex, fold and subx fields)
%
% The path looks like: footpath/save/fold/prex_subx.mat
%
% If there are redundant elements in both option and varargin, the following
% stategies will be used:
%    'svL': using save in option
%    'prex': prex1_prex2
%    'subx': subx1_subx2
%    'fold': fold1/fold2
%
% Example
%   input     -  option = {'svL', 0, 'prex', 'a1', 'fold', 'b1', 'subx', 'c1', 'path', ''};
%                varargins = {'svL', 1, 'redo', 'n', 'prex', 'a2', 'fold', 'b2', 'subx', 'c2'};
%   call      -  [svL, path, prex, fold, subx] = psSv(option, varargins{:});
%   output    -  svL = 0; 
%                path = 'footpath/b1/b2/a1_a2_c1_c2.mat';
%                prex = 'a1_a2';
%                fold = 'b1/b2';
%                subx = 'c1_c2';
%
% Input
%   option    -  1st option, 1 x (2 x m) (cell),
%                'svL' | 'prex' | 'fold' | 'subx' | 'path' | 'type'
%   varargin  -  2nd option
%                'svL' | 'prex' | 'fold' | 'subx'
%
% Output
%   svL      -  save level
%   path      -  path
%   prex      -  prefix
%   fold      -  fold
%   subx      -  subfix
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% parameter in option
svL1 = ps(option, 'svL', 0);
prex1 = ps(option, 'prex', '');
fold1 = ps(option, 'fold', '');
subx1 = ps(option, 'subx', '');
path  = ps(option, 'path', '');

% parameter in varargin
svL2 = ps(varargin, 'svL', 0);
prex2 = ps(varargin, 'prex', '');
fold2 = ps(varargin, 'fold', '');
subx2 = ps(varargin, 'subx', '');
type = ps(varargin, 'type', 'mat');

% integer in prex
if ~ischar(prex1)
    prex1 = sprintf('%d', prex1);
end
if ~ischar(prex2)
    prex2 = sprintf('%d', prex2);
end

svL = selectS(svL1, svL2, 'replace');

% not saved
if svL == 0
    path = [];
    prex = [];
    fold = [];
    subx = [];
    return;
end

global footpath; % specified in addPath.m
if ~isempty(path)
    prex = '';
    fold = '';
    subx = '';
    return;
end

% select one of the parameter
prex = selectS(prex1, prex2, 'join_', 'prex');
fold = selectS(fold1, fold2, 'join/', 'fold');

% if fold does not exist, create it
fold2 = sprintf('%s/save/%s', footpath, fold);
if ~isdir(fold2)
    mkdir(fold2);
end

% subx can be empty
if isempty(subx1) && isempty(subx2)
    subx = '';
    path = sprintf('%s/save/%s/%s.%s', footpath, fold, prex, type);
else
    subx = selectS(subx2, subx1, 'join_', 'subx');
    path = sprintf('%s/save/%s/%s_%s.%s', footpath, fold, prex, subx, type);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = selectS(s1, s2, type, info)
% Select the multiple string.%
%
% Input
%   s1      -  1st string
%   s2      -  2nd string
%   type    -  function type, 'replace' | 'join_' | 'joint/'
%   info    -  information of string, 'prex', 'fold', 'subx'
%
% Output
%   s       -  final string

flag1 = ~isempty(s1);
flag2 = ~isempty(s2);

if flag1 && flag2
    if strcmp(type, 'replace')
        s = s1;
    elseif strcmp(type, 'join_')
        s = [s1 '_' s2];
    elseif strcmp(type, 'join/')
        s = [s1 '/' s2];
    else
        error('unknown type');
    end

elseif flag1
    s = s1;

elseif flag2
    s = s2;

else
    error('Not enough parameters for ''%s''', info);
end
