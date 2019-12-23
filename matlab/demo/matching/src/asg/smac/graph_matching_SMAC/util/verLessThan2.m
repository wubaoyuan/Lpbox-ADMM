function tf=verLessThan2(strNumber);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

temp=version;
temp = getParts(temp);
strNumber = getParts(strNumber);

tf = (sign(temp - strNumber) * [1; .1; .01]) < 0;
% tf = (sign(toolboxParts - verParts) * [1; .1; .01]) < 0;

function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d')';
if length(parts) < 3
    parts(3) = 0; % zero-fills to 3 elements
end
