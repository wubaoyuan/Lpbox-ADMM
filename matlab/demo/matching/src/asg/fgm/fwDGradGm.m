function Gr = fwDGradGm(X)
% Compute the gradient of the original function part in Frank-Wolfe algorithm.
%
% Input
%   X       -  correspondence, n1 x n2
%           
% Output    
%   Gr      -  gradient, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-16-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variable
global KP KQ;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;
global G1s G2s H1s H2s;
global GXG HXH;
global GXGs HXHs;
global isGXG isHXH;

% ti = tic;
% if isGXG == 0
%     GXG = multGXH(IndG1T, X, IndG2);
%     isGXG = 1;
% end
% 
% if isHXH == 0
%     HXH = multGXH(IndH1T, X, IndH2);
%     isHXH = 1;
% end

% gradient: KP + H1 (G1' X G2 .* KQ) H2' + G1 (H1' X H2 .* KQ) G2'
% Gr = KP + multGXH(IndH1, GXG .* KQ, IndH2T) + multGXH(IndG1, HXH .* KQ, IndG2T);
% t1 = toc(ti);

%ti = tic;
%Xs = sparse(X);
if isGXG == 0
    GXGs = G1s' * X * G2s;
    isGXG = 1;
end
if isHXH == 0
    HXHs = H1s' * X * H2s;
    isHXH = 1;
end
Gr = KP + H1s * (GXGs .* KQ) * H2s' + G1s * (HXHs .* KQ) * G2s';
% t2 = toc(ti);
% fprintf('ratio: %f\n', t1 / t2);
% 
% equal('Gr', Gr, full(Grs));
