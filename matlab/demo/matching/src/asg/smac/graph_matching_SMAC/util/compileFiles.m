function filesError = compileFiles(filepaths);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

files = filepaths2files(filepaths);

dir_old = pwd;
n_old = length(dir_old);
for i=1:length(files)
    temp=files(i).filepath;
    if length(temp)>=n_old && strcmp(temp(1:n_old),dir_old)
        temp=temp(n_old+1:end);
    end
    files(i).filepathrel=temp;
end

mex_exts={'.cpp','.c'};
ind=find(strcmp({files.ext},mex_exts{1}) | strcmp({files.ext},mex_exts{2}));
files=files(ind);
isValid=isRealMexFiles({files.filepath});
files=files(isValid==1);

result=compileFiles_aux(files);

filesError=result([result.isError]==1);

nErrors=length(filesError);
if nErrors > 0
    disp([sprintf('\n'),'Error: There were ' num2str(nErrors) ' erroneous files during compilation']);
    for i=1:nErrors
        disp(filesError(i).file.filepathrel);
    end
else
    disp([sprintf('\n'),'Compilation of files succeded without error ']);%TODO: #files compiled
end
if nargout==0 && isempty(filesError)
    clear filesError;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=compileFiles_aux(files);
for i=1:length(files)
    result(i)=compileFile(files(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=compileFile(file);
[result.isErrorDuringCompilation,result.isCompilationFailed,result.message,result.compiledFilename]=compileFile_aux(file);
result.file=file;
result.isError=result.isErrorDuringCompilation || result.isCompilationFailed;
if result.isErrorDuringCompilation
    warning(['Error: uncaught error during compilation of ' file.filepathrel,' : ',result.message]);
else
    if result.isCompilationFailed
        disp('*********************************');
        disp(['Error: compilation of ' file.filepathrel,' failed',' : ',result.message]);
    else
        if ~isempty(result.compiledFilename)
            disp(['compiled ' file.filepathrel,' => ',result.compiledFilename]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isErrorDuringCompilation, isCompilationFailed, message, compiledFilename] = compileFile_aux(file)
oldDir = pwd;
isErrorDuringCompilation=0;
isCompilationFailed=0;
message='';
compiledFilename='';

try
    cd(file.path);
    
    if ispc
%         [isCompilationFailed,message]=mex_aux(file.name,[],[]);
        [isCompilationFailed,message]=mex_aux(file.name);
    else
%         [isCompilationFailed,message]=mex_aux(file.name,[],[]);
        [isCompilationFailed,message]=mex_aux(file.name);
    end
    if ~isCompilationFailed
        compiledFilename=changeExt(file.name,['.',mexext]);
    else
        isCompilationFailed=1;
    end
    cd(oldDir);
catch
    isErrorDuringCompilation=1;
    message=lasterr;    
    cd(oldDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isMexfiles = isRealMexFiles(filepaths)
for i = 1 : length(filepaths)
    isMexfiles(i) = isRealMexFile(filepaths{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isMexfile = isRealMexFile(filepath)
%TODO:do regexp to capture isolated word only
string = file2string(filepath);
isMexfile = any(strfind(string,'mexFunction'));

% sourceCode = textread(filepath,'%s','delimiter','\n','whitespace','');
% temp = strfind(sourceCode,'mexFunction');
% occurences = 0;
% for i=1:length(temp)
%     occurences = occurences + length(temp{i});
% end
% isMexfile=occurences>0;
