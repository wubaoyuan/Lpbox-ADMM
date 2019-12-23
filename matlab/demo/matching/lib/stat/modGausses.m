function [mes, Vars] = modGausses(k, varargin)
% Generate a set of Gaussian models.
%
% Input
%   k       -  #classes
%   varargin
%     dis   -  distribution of cluster mean, {'unif'} | 'gauss'
%
%     %% uniform distribution %%
%     meMi  -  minimum of mean, {[0; 0]} | d x 1
%     meMa  -  maximum of mean, {[10; 10]} | d x 1
%     rBd   -  bound of radius, {[.05, .15]} | 1 x 2
%
%     %% Gaussian distribution %%
%     me0   -  mean, dim x 1, {[5; 5]}
%     r0    -  radius, {1}
%     rBd   -  bound of radius, 1 x 2, {[.05, .15]}
%
%     maMiD -  flag of maximizing the minimum of pairwise distance, 'y' | {'n'}
%
% Output
%   mes     -  mean of Gaussian, 1 x k (cell), d x 1
%   Vars    -  variance of Gaussian, 1 x k (cell), d x d
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-19-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
dis = ps(varargin, 'dis', 'unif');
maMiD = ps(varargin, 'maMiD', 'n');

% generate mean
if strcmp(dis, 'unif')
    % hyper-parameter
    meMi = ps(varargin, 'meMi', [  0;   0]);
    meMa = ps(varargin, 'meMa', [ 10;  10]);
    rBd = ps(varargin, 'rBd', [.05, .13]);
    rMi = (meMa - meMi) * rBd(1);
    rMa = (meMa - meMi) * rBd(2);

    mes = smpUnif(meMi, meMa, k, 'maMiD', maMiD);
    
elseif strcmp(dis, 'gauss')
    % hyper-parameter
    me0 = ps(varargin, 'me0', [0; 0]);
    r0 = ps(varargin, 'r0', 1);
    Var0 = [r0, 0; 0, r0];
    rBd = ps(varargin, 'rBd', [.2, .5] * .8);
    rMi = [r0; r0] * rBd(1);
    rMa = [r0; r0] * rBd(2);

    mes = smpGauss(me0, Var0, k, 'maMiD', maMiD);

else
    error('unknown distribution: %s', dis);
end
mes = mdiv('horz', mes, ones(1, k));

% generate variance
rs = smpUnif(rMi, rMa, k);
as = pi * 2 * rand(1, k);
    
Vars = cell(1, k);
for c = 1 : k
    U = a2rot(as(c));
    Vars{c} = U' * diag(rs(:, c) .^ 2) * U;
end
