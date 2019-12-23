function  [result, timing] = eigs_optimized(A, B, k, options, Afun1, Afun2, Afun3)
% based on matlab's function eigs
% result = eigs_optimized_W_B(A,B,k,options);
%
% if isempty_B, A*X = lambda*X => OP = A and B = I
% if ~isempty_B, A*X = lambda*B*X is computed as :
% OP = R'\(A*(R\x)) where R is B's upper triangular Cholesky factor: B = R'*R.
% Finally, V = R\V returns the actual generalized eigenvectors of A and B.
%
% A is assumed symmetric, real, square
% B must be the same size as A and either a SP(D) matrix or its Cholesky factor
% (with or without permB)
%
% options : struct with fields :
% fast               :  [ 0 | 1  ] if 1, does not perform argument validity checking
% computeX           :  [ 0 | {1}] request or not to compute X
% warningConvergence :  [ 0 | {1}]
% sigma              :  [ {'LA'} | 'LM' ]
% n                  :  length of matrix represented by A if A is a function
% tol                :  default eps
% p                  :  default min(max(2*k,20),n)
% maxit              :  default max(300,ceil(2*n/max(p,1)))
% v0                 :  default zeros(n,1)
% cholB              :  [{0} | 1]
% permB              :  default 1:n
%
% output :
% result.X
% result.lambda
% result.flag
% result.nbA_times_X
% result.nbIterations
% result.nbEigenvaluesConverged
%
% timing.preprocessing
% timing.dsaupd
% timing.A_times_X
% timing.postprocessing
% timing.total
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

cputms = zeros(5, 1);
t0 = cputime; % start timing pre-processing

if isa(A, 'double')
    Amatrix = true;
    n = length(A);
else
    A = fcnchk(A);
    Amatrix = false;
    n = options.n;
    nbAfun = nargin - 4;
end
isempty_B = isempty(B);

fastMode = options.fastMode;
if fastMode
    rvec = int32(options.computeX);
    warningConvergence = options.warningConvergence;
    sigmaString = options.sigma;
    tol = options.tol;
    p = options.p;
    maxit = options.maxit;
    options.v0(1)=options.v0(1);%assures options.v0 will not be modified in place
    resid = options.v0;
    info = int32(1);
    if ~isempty_B
        cholB = options.cholB;
        if cholB
            permB = options.permB;
            RB = B;
            RBT = B';
        else
            if issparse(B)
                permB = symamd(B);
                [RB,pB] = chol(B(permB,permB));
            else
                [RB,pB] = chol(B);
            end
            if (pB == 0)
                RBT = RB';
            else
                error('BsizeMismatchAorNotSPDorNotChol');
            end
        end
    end

else
    if isfield(options,'sigma')
        if isequal(options.sigma,'LA') || isequal(options.sigma,'LM')
            sigmaString = options.sigma;
        else
            error('wrong argument');
        end
    else
        sigmaString = 'LA';
    end

    if isfield(options,'computeX') && options.computeX == 0
        rvec = int32(0);
    else
        rvec = int32(1);
    end
    if isfield(options,'warningConvergence') && options.warningConvergence == 0
        warningConvergence = 0;
    else
        warningConvergence = 1;
    end

    if ~options.issym
        error('issymA  == 0');
    end
    if (k>=n)
        error('k>=n')
    end

    if isfield(options,'tol')
        tol = options.tol;
    else
        tol = eps; % ARPACK's minimum tolerance is eps / 2 (DLAMCH's EPS)
    end
    if isfield(options,'p')
        p = options.p;
        if p>n || p<k
            error('p>n || p<k')
        end
    else
        p = min(max(2*k,20),n);
    end

    if isfield(options,'maxit')
        maxit = options.maxit;
    else
        maxit = max(300,ceil(2*n/max(p,1)));
    end

    if isfield(options,'v0')
        options.v0(1)=options.v0(1);%assures options.v0 will not be modified in place
        info = int32(1); % use resid as starting vector
        resid = options.v0;
    else
        info = int32(0); % use a random starting vector
        resid = zeros(n,1);
    end

    if ~isempty_B
        cholB = 0;
        if (isfield(options,'cholB') || isfield(options,'permB'))
            if isfield(options,'cholB')
                cholB = options.cholB;
                if isfield(options,'permB')
                    if issparse(B) && cholB
                        permB = options.permB;
                    else
                        warning('MATLAB:eigs:IgnoredOptionPermB', ...
                            ['Ignoring options.permB since B is not its sparse' ...
                            ' Cholesky factor.'])
                    end
                else
                    permB = 1:n;
                end
            end
        end

        if cholB
            RB = B;
            RBT = B';
        else
            if issparse(B)
                permB = symamd(B);
                [RB,pB] = chol(B(permB,permB));
            else
                [RB,pB] = chol(B);
            end
            if (pB == 0)
                RBT = RB';
            else
                error('BsizeMismatchAorNotSPDorNotChol');
            end
        end
    end

end

% Allocate outputs and ARPACK work variables
v = zeros(n,p);
% if isfield(options,'v')
%     v=options.v;    
%     assert(size(v,1)==n && size(v,2)==p);
% else
%     v = zeros(n,p);
% end

ldv = int32(size(v,1));
ipntr = int32(zeros(15,1));
workd = zeros(n,3);
lworkl = p*(p+8);
workl = zeros(lworkl,1);
lworkl = int32(lworkl);
d = zeros(k,1);

ido = int32(0); % reverse communication parameter
nev = int32(k); % number of eigenvalues requested
ncv = int32(p); % number of Lanczos vectors
iparam = int32(zeros(11,1));
iparam([1 3 7]) = int32([1 maxit 1]);
select = int32(zeros(p,1));
int32_n = int32(n);

cputms(1) = cputime - t0; % end timing pre-processing
nbA_times_X = 0;
while (ido ~= 99)

    t0 = cputime;

    arpackc('dsaupd', ido, 'I', int32_n, sigmaString , nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
%     timing_dsaupd(nbA_times_X+1)=cputime-t0;
    
    cputms(2) = cputms(2) + (cputime-t0);
    col1 = floor((double(ipntr(1))-1)/n) + 1;
    col2 = floor((double(ipntr(2))-1)/n) + 1;

    t0 = cputime; % start timing MATLAB OP(X)
    if Amatrix
        if isempty_B
            % OP = A*x
            workd(:,col2) = A * workd(:,col1);
        else
            % OP = R'\(A*(R\x))
            if issparse(B)
                workd(permB,col2) = RB \ workd(:,col1);
                workd(:,col2) = A * workd(:,col2);
                workd(:,col2) = RBT \ workd(permB,col2);
            else
                workd(:,col2) = RBT \ (A * (RB \ workd(:,col1)));
            end
        end
    else % A is not a matrix
        if isempty_B
            % OP = A*x
            switch nbAfun
                % workd(:,col2) = feval(A,workd(:,col1),varargin);
                case 0
                    workd(:,col2) = feval(A,workd(:,col1));
                case 1
                    workd(:,col2) = feval(A,workd(:,col1), Afun1);
                case 2
                    workd(:,col2) = feval(A,workd(:,col1), Afun1,Afun2);
                case 3
                    workd(:,col2) = feval(A,workd(:,col1), Afun1,Afun2,Afun3);
            end
        else
            % OP = R'\(A*(R\x))
            if issparse(B)
                workd(permB,col2) = RB \ workd(:,col1);
                % workd(:,col2) = feval(A,workd(:,col2),varargin);
                switch nbAfun
                    case 0
                        workd(:,col2) = feval(A,workd(:,col2));
                    case 1
                        workd(:,col2) = feval(A,workd(:,col2), Afun1);
                    case 2
                        workd(:,col2) = feval(A,workd(:,col2), Afun1,Afun2);
                    case 3
                        workd(:,col2) = feval(A,workd(:,col2), Afun1,Afun2,Afun3);
                end
                workd(:,col2) = RBT \ workd(permB,col2);
            else
                % workd(:,col2) = RBT \ feval(A,(RB\workd(:,col1)), varargin);
                switch nbAfun
                    case 0
                        workd(:,col2) = RBT \ feval(A,(RB\workd(:,col1)));
                    case 1
                        workd(:,col2) = RBT \ feval(A,(RB\workd(:,col1)), Afun1);
                    case 2
                        workd(:,col2) = RBT \ feval(A,(RB\workd(:,col1)), Afun1, Afun2);
                    case 3
                        workd(:,col2) = RBT \ feval(A,(RB\workd(:,col1)), Afun1, Afun2, Afun3);
                end
            end
        end
    end % if Amatrix
    cputms(3) = cputms(3) + (cputime-t0); % end timing MATLAB OP(X)
    nbA_times_X = nbA_times_X + 1;
    
    %{
    ki=p;
    sz=[66,100];
    X1=reshape(v(1:prod(sz),1:ki),[sz,ki]);    
    X1=reshape(v(1:prod(sz),1:ki),[sz,ki]);    
    if ~exist('X1old')
        X1old=X1;
    end
    Xd=cat(3,X1,X1-X1old);
    imagesc(showImages(Xd));
    drawnow;
    X1old=X1;
    0;
    %}
    
%     if nbA_times_X>maxit
%         break;
%     end
end % while (ido ~= 99)
nbIterations = double(ipntr(15));

t0 = cputime; % start timing post-processing

% result.v = v;
% result.v(1)=result.v(1);

arpackc( 'dseupd', rvec, 'A', select, d, v, ldv, 0, 'I', int32_n, sigmaString , nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info );
nconv = double(iparam(5));
if (nconv == 0)
    flag = 1;
    if warningConvergence
        warning('MATLAB:eigs:NoEigsConverged', 'None of the %d requested eigenvalues converged.',k)
    end
elseif (nconv < k)
    flag = 1;
    if warningConvergence
        warning('MATLAB:eigs:NotAllEigsConverged', 'Only %d of the %d requested eigenvalues converged.',nconv,k)
    end
else
    flag = 0;
end

indexesOrdered = [nconv:-1:1 , nconv+1:k];
if (rvec == 1)
    X = v(:,indexesOrdered);
    if ~isempty_B
        if issparse(B)
            X(permB,:) = RB \ X;
        else
            X = RB \ X;
        end
    end
end
lambda = d(indexesOrdered);
if nconv < k
    if ~isempty_B
        warning('Not yet supported. If you don''t care about the correct unconverged eigenvalues, you can ignore this');
    else
        for t = nconv+1:k,
            xt = v(:,t);
            if Amatrix
                yt = A * xt;
            else % A is not a matrix
                switch nbAfun
                    % workd(:,col2) = feval(A,workd(:,col1),varargin);
                    case 0
                        yt = feval(A,xt);
                    case 1
                        yt = feval(A,xt, Afun1);
                    case 2
                        yt = feval(A,xt, Afun1,Afun2);
                    case 3
                        yt = feval(A,xt, Afun1,Afun2,Afun3);
                end
            end % if Amatrix
            lambda(t) = xt'*yt;
        end
    end
end

if strcmp(options.sigma,'LA')
    [val,indexesOrdered] = sort(lambda,'descend');
elseif strcmp(options.sigma,'SA')
    [val,indexesOrdered] = sort(lambda,'ascend');
else
%     [val,indexesOrdered] = sort(lambda,'descend');
    error('TODO');
end
lambda = lambda(indexesOrdered);
if (rvec == 1)
    X = X(:,indexesOrdered);
end

% lambda = d(k:-1:1);
% if (rvec == 1)
%     X = v(:,k:-1:1);
%     if ~isempty_B
%         if issparse(B)
%             X(permB,:) = RB \ X;
%         else
%             X = RB \ X;
%         end
%     end
% end


cputms(4) = cputime-t0; % end timing post-processing
cputms(5) = sum(cputms(1:4)); % total time

result.lambda = lambda;
if (rvec == 1)
    result.X = X;
end

result.flag = flag;
result.nbEigenvaluesConverged = nconv;

timing.nbA_times_X = nbA_times_X;
timing.nbIterations = nbIterations;

timing.preprocessing = cputms(1);
timing.dsaupd = cputms(2);
timing.A_times_X = cputms(3);
timing.postprocessing = cputms(4);
timing.total = cputms(5);



