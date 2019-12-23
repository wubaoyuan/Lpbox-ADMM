function [options, nbEigenvectors] = getDefaultOptionsEigs(n, nbEigenvectors)

options.issym = 1;
options.maxit = 500;
options.tol = 1e-3;
options.v0 = ones(n, 1);
options.p = nbEigenvectors * 2;

if nbEigenvectors == 1
    options.p = 3;
end
options.p = min(round(options.p), n);

options.fastMode = 1;
options.computeX = 1;
options.warningConvergence = 1;
options.sigma = 'LA';
options.n = n;
nbEigenvectors = min(nbEigenvectors, options.p - 1);
