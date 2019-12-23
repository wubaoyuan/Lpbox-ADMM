function [Ctx, Pt] = shpCtx(Pt0, par)
% Compute shape context from a set of points.
%
% References
%   S. Belongie, J. Malik and J. Puzicha
%   "Shape Matching and Object Recognition Using Shape Contexts", PAMI, 2002
%
% Input
%   Pt0      -  original point set, 2 x nPt0
%   par      -  function parameter
%     smp    -  flag of re-sampling, {'y'} | 'n'
%     nPt    -  #points after re-sampling, {100}
%     nBinT  -  #bins in circle, {12}
%     nBinR  -  #bins in radius, {5}
%
% Output
%   Ctx      -  shape context, (nBinT x nBinR) x nPt
%   Pt       -  new point set after sampling, 2 x nPt
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-26-2010
%   modify   -  Feng Zhou (zhfe99@gmail.com), 12-09-2012

% function parameter
isSmp = psY(par, 'smp', 'y');
nPt = ps(par, 'nPt', 100);
nBinT = ps(par, 'nBinT', 12);
nBinR = ps(par, 'nBinR', 5);

% dimension
prIn('shpCtx', 'isSmp %d, nPt %d, nBinT %d, nBinR %d', ...
     isSmp, nPt, nBinT, nBinR);

% re-sampling
if isSmp
    Pt = ctrSmp(Pt0, nPt);
else
    Pt = Pt0;
    nPt = size(Pt, 2);
end

% shape context
r_inner = 1 / 8;
r_outer = 2;
Ctx = sc_compute(Pt, zeros(1, nPt), [], nBinT, nBinR, r_inner, r_outer, zeros(1, nPt));

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pt = ctrSmp(Pt0, nPt)
% Jitendra's method for sampling points.
%
% Input
%   Pt0  -  position of original points, 2 x nPt0
%   nPt  -  #sampling points
%
% Output
%   Pt   -  position of new points, 2 x nPt

% dimension
nPt0 = size(Pt0, 2);
k = 3;
nPt1 = min(k * nPt, nPt0);

% randomly re-order
ind0 = randperm(nPt0);
ind0 = ind0(1 : nPt1);
Pt1 = Pt0(:, ind0);

% pairwise distance
D1 = conDst(Pt1, Pt1);
D1 = D1 + diag(Inf * ones(nPt1, 1));

% remove points that are too close
while nPt1 ~= nPt
   % find closest pair
   [a, b] = min(D1);
   [c, d] = min(a);
   I = b(d);
   J = d;

   % remove one point
   Pt1(:, J) = [];
   D1(:, J) = [];
   D1(J, :) = [];
   nPt1 = nPt1 - 1;
end

Pt = Pt1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BH = sc_compute(Bsamp, Tsamp, mean_dist, nbins_theta, nbins_r, r_inner, r_outer, out_vec)
% Compute shape histogram.
%
% Input
%   Bsamp      -  position of points, 2 x nsamp
%   Tsamp      -  1 x nsamp (tangent theta)
%   mean_dist  -  mean distance, used for length normalization
%   out_vec    -  1 x nsamp (0 for inlier, 1 for outlier)
%                 outliers are not counted in the histograms, but they do get assigned a histogram
%
% Output
%  BH          -  histogram, nbins x nsamp

% dimension
nsamp = size(Bsamp, 2);
in_vec = out_vec == 0;

% compute r, theta arrays
D = conDst(Bsamp, Bsamp);
r_array = real(sqrt(D));
theta_array_abs = atan2(Bsamp(2, :)' * ones(1, nsamp) - ones(nsamp, 1) * Bsamp(2, :), ...
                        Bsamp(1, :)' * ones(1, nsamp) - ones(nsamp, 1) * Bsamp(1, :))';
theta_array = theta_array_abs - Tsamp' * ones(1, nsamp);

% create joint (r,theta) histogram by binning r_array and theta_array

% normalize distance by mean, ignoring outliers
if isempty(mean_dist)
   tmp = r_array(in_vec,:);
   tmp = tmp(:, in_vec);
   mean_dist = mean(tmp(:));
end
r_array_n = r_array / mean_dist;

% use a log scale for binning the distances
r_bin_edges = logspace(log10(r_inner), log10(r_outer), 5);
r_array_q = zeros(nsamp, nsamp);
for m = 1 : nbins_r
   r_array_q = r_array_q + (r_array_n < r_bin_edges(m));
end

% flag all points inside outer boundary
fz = r_array_q > 0;

% put all angles in [0, 2pi) range
theta_array_2 = rem(rem(theta_array, 2 * pi) + 2 * pi, 2 * pi);

% quantize to a fixed set of angles (bin edges lie on 0, (2 * pi) / k, ..., 2 * pi
theta_array_q = 1 + floor(theta_array_2 / (2 * pi / nbins_theta));

nbins = nbins_theta * nbins_r;
BH = zeros(nbins, nsamp);
for n = 1 : nsamp
    fzn = fz(n, :) & in_vec;
    Sn = sparse(theta_array_q(n, fzn), r_array_q(n, fzn), 1, nbins_theta, nbins_r);
    BH(:, n) = Sn(:);
end
