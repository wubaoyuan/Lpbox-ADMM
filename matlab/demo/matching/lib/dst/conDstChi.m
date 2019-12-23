function D = conDstChi(Hist1, Hist2, varargin)
% Calculate the chi-square distance of two histograms.
%   \chi(P, Q) = .5 * \sum_i (P_i - Q_i)^2 / (P_i + Q_i)
%
% Remark
%   Notice that chi-square distance is not a metric (i.e., not satisfied triangle inequality).
%
% Input
%   Hist1   -  1st histogram, nB x n1
%   Hist2   -  2nd histogram, nB x n2
%   varargin
%     norm  -  normalization flag, {'y'} | 'n'
%
% Output
%   D       -  squared distance matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-11-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-07-2012

% function option
isNorm = psY(varargin, 'norm', 'y');

% dimension
[nB, n1] = size(Hist1);
n2 = size(Hist2, 2);

% normalization
if isNorm
    hist1 = sum(Hist1) + eps;
    Hist1 = Hist1 ./ repmat(hist1, nB, 1);

    hist2 = sum(Hist2) + eps;
    Hist2 = Hist2 ./ repmat(hist2, nB, 1);
end

D = zeros(n1, n2);
for iB = 1 : nB
    H1 = repmat(Hist1(iB, :)', 1, n2);
    H2 = repmat(Hist2(iB, :), n1, 1);
    D = D + (H1 - H2) .^ 2 ./ (H1 + H2 + eps);
end
D = D * .5;

% squared
%D = D .^ 2;
