function [c_tps, d_tps, w] = rpmTranOpt(algT, lamda1, lamda2, sigma, x, vy)
% Calculate transformation for RPM.
% 
% Input
%   algT    -  transformation type, 'tps' | 'rbf'
%   lamda1  -
%   lamda2  -
%   sigma   -
%   x       -
%   vy      -
%
% Output
%   ** for TPS **
%   c_tps   -  TPS weight, n x 3
%   d_tps   -  affine transformation, 3 x 3
%   ** for RBF **
%   w       -  weight, n x 3
%
% History
%   create  -  Anand Rangarajan (anand@noodle.med.yale.edu), 04-27-2000
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-20-2012

if strcmp(algT, 'tps')
    [c_tps, d_tps] = rpmTranOptTps(x, vy, lamda1, lamda2);
    w = [];

elseif strcmp(algT, 'rbf')
    c_tps = [];
    d_tps = [];
    [~, w] = crbf_gen(x, vy, x, lamda1, lamda2, sigma);

else
    error('unknown transformation type: %s', algT);
end
