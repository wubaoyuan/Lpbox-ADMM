function objGm = pathDObjGm(X)
% Compute objective for the path-following algorithm.
%
% Input
%   X         -  correspondence, n1 x n2
%             
% Output      
%   objGm     -  J_gm
%
% History     
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-26-2012
%   modify    -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global
global KP KQ;
global HH1 HH2;
global G1s G2s H1s H2s HH1s HH2s;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T IndHH1 IndHH2;
global indG1 indG2 indH1 indH2;
global QQ1 QQ2 gamma;

%ti = tic;
% compute: J_gm = trace(KP' X) + trace(KQ' (G1' X G2 .* H1' X H2))
% tmp1 = multTr(KP .* X);
% tmp2 = multTr(KQ .* multGXH(IndG1T, X, IndG2) .* multGXH(IndH1T, X, IndH2));
% objGm = tmp1 + tmp2;

% comptue: J_con = trace((Q1' Q1)' (G1' X X' G1 .* H1' H1)) 
%                + trace((Q2' Q2)' (G2' X X' G2 .* H2' H2))
% XX = X * X';
% tmp1 = multTr(QQ1, multGXH(IndG1T, XX, IndG1), HH1);
% tmp2 = multTr(QQ2, multGXH(IndG2T, XX, IndG2), HH2);
% objCon = tmp1 + tmp2;

% objVex = objGm - .5 * objCon;
% objCav = objGm + .5 * objCon;

% linear interoplation
% obj = objGm + (alp - .5) * objCon;
% obj2 = (1 - alp) * objVex + alp * objCav;
% if ~equal('obj', obj, obj2, 'pr', 'n')
%     pr('problem in computing objective');
% end
% t1 = toc(ti);

% ti = tic;
Xs = sparse(X);
GXGs = G1s' * Xs * G2s;
HXHs = H1s' * Xs * H2s;

tmp1 = sum(sum(KP .* Xs));
tmp2 = sum(sum(GXGs .* HXHs .* KQ));
objGm = tmp1 + tmp2;
