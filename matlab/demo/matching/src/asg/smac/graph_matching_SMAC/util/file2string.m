function string = file2string(file);
%TODO: deal with differences ('\r\n' vs \n vs \r etc)
%{
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

string = file2string(which('file2strings'));
%}


string = file2string_fscanf(file);

% strings=file2strings(file);
% string2=strings2string(strings);
% isequal(string,string2)
0;


function string = file2string_fscanf(file);
fid = fopen(file);
assert(fid~=-1);
string = fscanf(fid, '%c', inf);
fclose(fid);

% string(string==sprintf('\r'))=[];
% doesn't work when isolated \r

string2=[string(2:end),'a'];
string(string==sprintf('\r') & string2==sprintf('\n'))=[];
string(string==sprintf('\r'))=sprintf('\n');

0;
% string=regexprep(string,'\r',''); %slow


% A = fread(fid)
% string = textscan(fid,'%s');
% string = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
% fclose(fid);



%{
bufsize = 4095;
if ~exist(file)
    error('file does not exist');
end
while 1
    try
        string = textread(file, '%s','bufsize',bufsize);
        string = string{1};
        break;
    catch
        bufsize = 2*bufsize;
        if bufsize > 16*4095
            error('bufsize > 16*4095');
        end
    end
end
%}
