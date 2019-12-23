function files = dir2(dirInput);
% Timothee Cour, 21-Apr-2008 17:31:23

files = dirRec(dirInput);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function files = dirRec(dirInput)
index = 0;
dirRoot = dirInput;
files = [];
dirRec_aux(0, dirInput, '');

    function dirRec_aux(parent,directory,directoryRel);
        pause(0);

        files0 = dir(directory);
        files0 = files0(3:end);

        for i=1:length(files0)
            index=index+1;
            file_i=newFile(files0(i),parent,directory,directoryRel);
            if isempty(files)
                files=file_i;
            else                
                files(index)=file_i;
            end
            if file_i.isdir
                dir_i=fullfile(directory,file_i.name);
                dir_rel_i=fullfile(directoryRel,file_i.name);
                dirRec_aux(index,dir_i,dir_rel_i);
            end

        end

    end

    function file=newFile(file,parent,dirParent,dirParentRel);
        [ignore,name0,ext]=fileparts(file.name);
        file.path=dirParent;
        file.name0=name0;
        file.ext=ext;
        file.filepath=fullfile(dirParent,file.name);
        file.filepathrel=fullfile(dirParentRel,file.name);
        file.parent=parent;
        file.index=index;
    end
end
