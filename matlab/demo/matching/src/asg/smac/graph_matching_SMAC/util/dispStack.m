function dispStack(stack);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<1
    if verLessThan2('7.1')
        disp('no stack available in this version');
        return;
    end
    temp=lasterror;
    disp(['handled error: ',temp.message]);
    stack=temp.stack;
end
% disp('stack: ');
for i=1:length(stack)
    disp(['In ','<a href="matlab: opentoline(''',stack(i).file,''',',num2str(stack(i).line),');">',stack(i).name,' at ',num2str(stack(i).line),'</a>']);    
end

% com.mathworks.mlservices.MLEditorServices.toString(stack(i).file);
% com.mathworks.mlservices.MLEditorServices.builtinGetNumOpenDocuments()

