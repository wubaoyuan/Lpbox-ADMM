function y = std(varargin)
%STD Standard deviation.
%   For vectors, Y = STD(X) returns the standard deviation.  For matrices,
%   Y is a row vector containing the standard deviation of each column.  For
%   N-D arrays, STD operates along the first non-singleton dimension of X.
%
%   STD normalizes Y by (N-1), where N is the sample size.  This is the
%   sqrt of an unbiased estimator of the variance of the population from
%   which X is drawn, as long as X consists of independent, identically
%   distributed samples.
%
%   Y = STD(X,1) normalizes by N and produces the square root of the second
%   moment of the sample about its mean.  STD(X,0) is the same as STD(X).
%
%   Y = STD(X,FLAG,DIM) takes the standard deviation along the dimension
%   DIM of X.  Pass in FLAG==0 to use the default normalization by N-1, or
%   1 to use N.
%
%   Example: If X = [4 -2 1
%                    9  5 7]
%     then std(X,0,1) is [3.5355 4.9497 4.2426] and std(X,0,2) is [3.0
%                                                                  2.0]
%   Class support for input X:
%      float: double, single
%
%   See also COV, MEAN, VAR, MEDIAN, CORRCOEF.
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.25.4.2 $  $Date: 2004/03/09 16:16:30 $

% Call var(x,flag,dim) with as many of those args as are present.
y = sqrt(var2(varargin{:}));
