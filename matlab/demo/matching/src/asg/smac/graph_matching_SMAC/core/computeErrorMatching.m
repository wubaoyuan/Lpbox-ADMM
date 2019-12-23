function errorMatching = computeErrorMatching(Xd, target, E12);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if isempty(target)
    errorMatching=[];
    return;
end

%TODO:handle when sum(target)~=n1
[n1, n2] = size(E12);
nMin = min(n1, n2);

errorMatching.nbErrors = sum(vec(Xd ~= target)) / 2; % because don't count an error twice in case of permutation
errorMatching.errorRate = errorMatching.nbErrors / nMin;
errorMatching.randomErrorRate = (nMin - 1) / nMin;
