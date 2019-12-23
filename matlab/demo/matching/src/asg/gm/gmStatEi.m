function [xEs, hstEs] = gmStatEi(eis)
% Compute the statistics of the kernel used in graph matching.
%
% Input
%
% Output
% 
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-25-2013

% dimension
len = length(eis);

% eigen-value of K
vEs = zeros(1, len);
for i = 1 : len
    ei = eis{i};
    
    if ~isreal(ei)
        pr('%d, complex eigen-value', i);
        ei = real(ei);
    end
    
    mi = min(ei);
    ma = max(ei);
    vEs(i) = mi / (ma + eps);
end

nBin = 100;
stE = linspace(-1, 0, nBin + 1);
hstEs = x2hst(vEs, stE);
xEs = stE(1 : end - 1) + 1 / nBin;
