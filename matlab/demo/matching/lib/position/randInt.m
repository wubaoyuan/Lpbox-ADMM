function i = randInt(n)
% Randomly pick one integer in (1 : n) with the given probability.
%
% Input
%   n       -  number of sampling times
%
% Output
%   i    -  integer
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-04-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 06-10-2012

i = round(rand(1) * (n - 1)) + 1;