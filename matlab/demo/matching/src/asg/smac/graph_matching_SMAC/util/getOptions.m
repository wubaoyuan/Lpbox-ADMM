function options2=getOptions(optionsDefault,options);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

options2=optionsDefault;
if nargin<2 || isempty(options)
    return;
end
names = fieldnames(options);
for i=1:length(names)
    options2.(names{i})=options.(names{i});
end
