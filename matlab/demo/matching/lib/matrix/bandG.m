function sigma = bandG(D, nei)
% Compute the bandwidth of Gauss kernel.
%
% Input
%   D       -  squared distance matrix, n1 x n2
%              If D == [], then a distance matrix named "DG" has been put as a global variable.
%   nei     -  #nearest neighbour to compute the kernel bandwidth, {.1}
%              NaN: sigma == 1 which indicates nomalization doesn't apply
%              0:   sigma == 0 which indicates the binary kernel
%
% Output
%   sigma   -  kernel bandwidth (variance)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% global
global DG;
isDG = isempty(D);

if isnan(nei)
    sigma = 1;
    return;
end

if nei == 0
    sigma = 0;
    return;
end

% dimension
if isDG
    n = size(DG, 1);
else
    n = size(D, 1);
end
m = min(max(1, floor(n * nei)), n);

% nearest neighbours
if isDG
    Dsorted = sort(DG);
else
    Dsorted = sort(D);
end

D2 = real(sqrt(Dsorted(1 : m, :)));
sigma = sum(D2(:)) / (n * m);
