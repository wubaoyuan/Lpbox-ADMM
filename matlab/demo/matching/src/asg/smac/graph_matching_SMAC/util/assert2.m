function assert2(condition, msgid, str);
%{
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

problem with builtin assert: assert(0) not defined
%}

if ~condition
    dispStack();
    error(['assert() was false']);
%     error(['assert(',inputname,') is false']);
end

%{
function assert2(condition);
if verLessThan2('7.4')
    condition=varargin{1};
    if ~condition
        error(['assert() was false']);
        error(['assert(',inputname,') is false']);
    end
else
    varargin{1}=logical(varargin{1});
    [varargout{1:nargout}] = builtin('assert', varargin{:});
end
%}

