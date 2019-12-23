function val=computeObjectiveValue(W,x,isNormalized);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

x=x(:);
if ~isequal(W,tril(W))
    W=tril(W);
end

val=x'*mex_w_times_x_symmetric_tril(x,W);

if nargin<3
    isNormalized=0;
end

if isNormalized
    val=val/(x'*x);
end

