function KQ = knlEgNor(KQ0, par)
% Normalize edge affinity.
%
% Reference
%   T. Cour, P. Srinivasan, and J. Shi, "Balanced graph matching", In NIPS, 2006
%
% Input
%   KQ0      -  original edge affinity, m1 x m2
%   par      -  parameter
%     nor    -  flag of whether to normalize edge affinity, 'y' | {'n'}
%     nItMa  -  #maximum iteration, {10}
%     th     -  threshold, {1e-7}
%
% Output
%   KQ       -  new edge affinity, m1 x m2
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-13-2011

% function parameter
isNor = psY(par, 'nor', 'n');
nItMa = ps(par, 'nItMa', 10);
th = ps(par, 'th', 1e-7);

if ~isNor
    KQ = KQ0;
    return;
end

% dimension
[m1, m2] = size(KQ0);

tolC = 1e-3;
KQ = bistocNormalize_slack(KQ0, tolC);
KQ = KQ / max(max(KQ));
