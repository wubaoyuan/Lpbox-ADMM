function disp2(expr,result);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin>=2
    disp_aux(expr,result);
    return;
end

try
    % result = evalin('base',expr);
    result = evalin('caller',expr);
catch
    assert(0);
    return;
    %     result=expr;
    % %     inputname(1)
    %     expr=inputname(1);
    %     result = evalin('caller',expr);
end
disp_aux(expr,result);


function disp_aux(expr,result);
fprintf([expr,' : ']);
if isvector(result) % && size(result,1)==1
    result=result(:)';
%     result
    %     fprintf([mat2str(result,5),'\n']);
    fprintf([num2str(result),'\n']);
else
    fprintf('\n');
    disp(result);
end
