function [cs, k, dict, dict0] = classDict(c0s)
% Build a dictionary which maps the class id into the real underlying class.
% Given a set of labels, we usually use a part of them.
% Therefore, the built CMap matrix maps the used class to real underlying class.
%
% Example
%   input    -  c0s = [2 2 3];
%   call     -  [cs, k, dict, dict0] = classDict(c0s);
%   output   -  cs = [1 1 2];
%               k = 2;
%               dict = [2 3];
%               dict0 = [0 1 2];
%
% Input
%   c0s      -  original class labels, 1 x m
%
% Output
%   cs       -  new class labels, 1 x m
%   k        -  #different classes
%   dict     -  dictionary (new -> original), 1 x k
%   dict0    -  dictionary (original -> new), 1 x cMa
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
k = 0;
m = length(c0s);
cMa = max(c0s);

% original -> new
dict0 = zeros(1, cMa);
for i = 1 : m
    c0 = c0s(i);
    if dict0(c0) == 0
        k = k + 1;
        dict0(c0) = k;
    end
end
cs = dict0(c0s);

% new -> original
dict = zeros(1, k);
for i = 1 : m
    dict(cs(i)) = c0s(i);
end
