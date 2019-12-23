function vx = rpmWarp(algT, x, c_tps, d_tps, w, sigma)
% Apply transformation for RPM.
%
% Input
%   algT
%   x
%   c_tps
%   d_tps
%   w
%   sigma
%
% Output
%   vx
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-20-2012

if strcmp(algT, 'tps')
    vx = ctps_warp_pts(x, x, c_tps, d_tps);

elseif strcmp(algT, 'rbf')
    vx = crbf_warp_pts(x, x, w, sigma);

else
    error('unknown transformation type: %s', algT);
end
