function [errorCode,message] = mex_silent(isVerbose,varargin)
%MEX Compile MEX-function. 
%   MEX [option1 ... optionN] sourcefile1 [... sourcefileN]
%       [objectfile1 ... objectfileN] [libraryfile1 ... libraryfileN]
%   
%   Description:
%     MEX compiles and links source files into a shared library called a
%     MEX-file, executable from within MATLAB. The resulting file has a
%     platform-dependent extension, as shown in the table below:
%   
%       solaris         - .mexsol
%       hpux            - .mexhpux
%       glnx86          - .mexglx
%       glnxi64         - .mexi64
%       Mac OS X        - .mexmac
%       Windows         - .dll
%   
%     The first file name given (less any file name extension) will be the name
%     of the resulting MEX-file. Additional source, object, or library files can
%     be given to satisfy external references. On Windows, either C or Fortran,
%     but not both, may be specified. On UNIX, both C and Fortran source files
%     can be specified when building a MEX-file. If C and Fortran are mixed, the
%     first source file given determines the entry point exported from the
%     MEX-file (MATLAB loads and runs a different entry point symbol for C or
%     Fortran MEX-files).
%   
%     Both an options file and command line options affect the behavior of MEX.
%     The options file contains a list of variables that are passed as arguments
%     to various tools such as the compiler, linker, and other platform-
%     dependent tools (such as the resource linker on Windows). Command line
%     options to MEX may also affect what arguments are passed to these tools,
%     or may control other aspects of MEX's behavior.
%   
%   Command Line Options:
%     Options available on all platforms:
%   
%     -ada <sfcn.ads>
%         Use this option to compile a Simulink S-function written in Ada, where
%         <sfcn.ads> is the Package Specification for the S-function. When this
%         option is specified, only the -v (verbose), -g (debug) and 
%         -I<pathname> options are relevant. All other options are ignored. See
%         $MATLAB/simulink/ada/examples/README for examples and info on
%         supported compilers and other requirements.
%     -argcheck
%         Add argument checking. This adds code so that arguments passed
%         incorrectly to MATLAB API functions will cause assertion failures.
%         Adds -DARGCHECK to the C compiler flags, and adds
%         $MATLAB/extern/src/mwdebug.c to the list of source files. (C functions
%         only).
%     -c
%         Compile only. Creates an object file but not a MEX-file.
%     -D<name>
%         Define a symbol name to the C preprocessor. Equivalent to a
%         "#define <name>" directive in the source.
%     -D<name>#<value>
%         Define a symbol name and value to the C preprocessor. Equivalent to a
%         "#define <name> <value>" directive in the source.
%     -f <optionsfile>
%         Specify location and name of options file to use. Overrides MEX's
%         default options file search mechanism.
%     -g
%         Create a debuggable MEX-file. If this option is specified, MEX appends
%         the value of options file variables ending in DEBUGFLAGS with their
%         corresponding base variable. (For example, the value of LINKDEBUGFLAGS
%         would be appended to the LINKFLAGS variable before calling the
%         linker.) This option also disables MEX's default behavior of
%         optimizing built object code.
%     -h[elp]
%         Print this message.
%     -I<pathname>
%         Add <pathname> to the list of directories to search for #include
%         files.
%     -inline
%         Inline matrix accessor functions (mx*). The MEX-function generated
%         may not be compatible with future versions of MATLAB.
%     -n
%         No execute mode. Print out any commands that MEX would otherwise have
%         executed, but do not actually execute any of them.
%     -O
%         Optimize the object code by including the optimization flags listed in
%         the options file. If this option is specified, MEX appends the value
%         of options file variables ending in OPTIMFLAGS with their
%         corresponding base variable. (For example, the value of LINKOPTIMFLAGS
%         would be appended to the LINKFLAGS variable before calling the
%         linker.) Note that optimizations are enabled by default, are disabled
%         by the -g option, but are reenabled by -O.
%     -outdir <dirname>
%         Place all output files in directory <dirname>.
%     -output <resultname>
%         Create MEX-file named <resultname> (an appropriate MEX-file extension
%         is automatically appended). Overrides MEX's default MEX-file naming
%         mechanism.
%     -setup
%         Interactively specify the compiler options file to use as default for
%         future invocations of MEX by placing it in the user profile
%         directory returned by PREFDIR. When this option is specified, no 
%         other command line input is accepted.
%     -U<name>
%         Remove any initial definition of the C preprocessor symbol <name>.
%         (Inverse of the -D option.)
%     -v
%         Print the values for important internal variables after the options
%         file is processed and all command line arguments are considered.
%         Prints each compile step and final link step fully evaluated to see
%         which options and files were used. Very useful for debugging.
%     -V5
%         Compile a MATLAB version 5-style MEX-file. This option is intended as
%         an aid to migration, and is not recommended as a permanent solution.
%     <name>#<value>
%         Override an options file variable for variable <name>. See the
%         platform-dependent discussion of options files below for more details.
%         This option is processed after the options file is processed and all
%         command line arguments are considered.
%   
%   Additional options available on Windows platforms:
%   
%     @<rspfile>
%         Include contents of the text file <rspfile> as command line arguments
%         to MEX.
%   
%   Additional options available on Unix platforms:
%   
%     -<arch>
%         Assume local host has architecture <arch>. Possible values for <arch>
%         include sol2, hpux, hp700, alpha, ibm_rs, sgi, and glnx86.
%     -D<name>=<value>
%         Define a symbol name and value to the C preprocessor. Equivalent to a
%         "#define <name> <value>" directive in the source.
%     -fortran
%         Specify that the gateway routine is in Fortran. This will override
%         what the script normally assumes, which is that the first source file
%         in the list is the gateway routine.
%     -l<name>
%         Link with object library "lib<name>" (for "ld(1)").
%     -L<directory>
%         Add <directory> to the list of directories containing object-library
%         routines (for linking using "ld(1)").
%     <name>=<value>
%         Override an options file variable for variable <name>. See the
%         platform-dependent discussion of options files below for more details.
%   
%   Options File Details:
%     On Windows:
%       The options file is written as a DOS batch file. If the -f option is not
%       used to specify the options file name and location, then MEX searches
%       for an options file named mexopts.bat in the following directories: the
%       current directory, then the user profile directory (returned by the
%       PREFDIR function), and lastly the directory specified by [matlabroot
%       '\bin\win32\mexopts']. Any variable specified in the options file
%       can be overridden at the command line by use of the <name>#<value>
%       command line argument. If <value> has spaces in it, then it should be
%       wrapped in double quotes (e.g., COMPFLAGS#"opt1 opt2"). The definition
%       can rely on other variables defined in the options file; in this case
%       the variable referenced should have a prepended "$" (e.g.,
%       COMPFLAGS#"$COMPFLAGS opt2").
%   
%       Note: The options files in $MATLAB\bin\mexopts named *engmatopts.bat are
%       special case options files that can be used with MEX (via the -f option)
%       to generate stand-alone MATLAB Engine and MATLAB MAT-API executables.
%       Such executables are given a ".exe" extension.
%   
%     On UNIX:
%       The options file is written as a UNIX shell script. If the -f option 
%       is not used to specify the options file name and location, then MEX 
%       searches for an options file named mexopts.sh in the following 
%       directories: the current directory, then the user profile directory 
%       (returned by the PREFDIR function), and lastly the directory specified 
%       by [matlabroot '/bin']. Any variable specified in the options file can
%       be overridden at the command line by use of the <name>=<def> command 
%       line argument. If <def> has spaces in it, then it should be wrapped in 
%       single quotes (e.g., CFLAGS='opt1 opt2'). The definition can rely on 
%       other variables defined in the options file; in this case the variable 
%       referenced should have a prepended "$" (e.g., CFLAGS='$CFLAGS opt2').
%   
%       Note: The options files in $MATLAB/bin named engopts.sh and matopts.sh
%       are special case options files that can be used with MEX (via the -f
%       option) to generate stand-alone MATLAB Engine and MATLAB MAT-API
%       executables. Such executables are not given any default extension.
%   
%   Examples:
%       The following command will compile "myprog.c" into "myprog.mexsol" (when
%       run under Solaris):
%   
%         mex myprog.c
%   
%       When debugging, it is often useful to use "verbose" mode as well
%       as include symbolic debugging information:
%   
%         mex -v -g myprog.c
%
%   See also MEXDEBUG, JAVA, PCODE, PERL, PREFDIR
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.102.4.3 $

% Normal behavior for MEX is to return an error code only if one is
% requested (nargout > 0).  If, instead, it is better to have it 
% return an error code regardless of nargout, then set the global
% variable MEX_RETURN_RESULT_IN_ANS to 1.  MEX will then return its
% result such that MATLAB puts that error code into the "ans" 
% default variable.
global MEX_RETURN_RESULT_IN_ANS
c = computer;

message = '';
if isunix

    [mexname, setup] = get_mex_opts(varargin{:});
    
    if ~isempty(mexname)
        [loaded_m, loaded_mex] = inmem;
    
        if ~isempty(loaded_mex)
            clear_mex_file(mexname);
        end

    end
    
    args = [' "ARCH=' lower(computer) '"'];
	if (nargin > 0)
		args = [args sprintf(' "%s"', varargin{:})];
	end
% 	errCode = unix([matlabroot '/bin/mex' args]);
	[errCode,message] = unix([matlabroot '/bin/mex' args]);
    if isVerbose
        disp(message);
    end
elseif strncmp(c, 'PC', 2)
    
    mexname = get_mex_opts(varargin{:});
    matlab_bin_location=[matlabroot '\bin'];

    %% changed for matlab 7.1
    if exist([matlab_bin_location '\win32\mex.bat'],'file') 
        if exist([matlab_bin_location '\win32\mex.pl'],'file')
            matlab_bin_location=[matlabroot '\bin\win32'];
        end
    end
    %{
    if exist([matlab_bin_location '\win32\mex.bat'],'file') 
        matlab_bin_location=[matlabroot '\bin\win32'];
    end
    %}
    
    if ~isempty(mexname)
        [loaded_m, loaded_mex] = inmem;
        if ~isempty(loaded_mex)
            clear_mex_file(mexname);
        end
    end
    
    % Loop over all the arguments. Put extra quotes around any that
    % contain spaces.
    
    for i=1:prod(size(varargin))
        if (find(varargin{i} == ' '))
            varargin{i} = [ '"' varargin{i} '"' ];
        end
    end
    
    % Format the mex command
    cmdargs = ['-called_from_matlab -matlab "' matlabroot '" ' sprintf(' %s', varargin{:})];
    if (any(matlab_bin_location == ' '))
        quote_str = '"';
    else
        quote_str = '';
    end

    cmdtool = [quote_str matlabroot '\sys\perl\win32\bin\perl.exe' quote_str ' ' ... 
               quote_str matlab_bin_location '\mex.pl' quote_str];
    [cmd, rspfile] = make_rsp_file(cmdtool, cmdargs);
    try
        [errCode,message] = dos(cmd);
        if isVerbose
            disp(message);
        end
    catch
        if isVerbose
            disp(lasterr);
        end
        errCode = 1; % failure
    end
    delete(rspfile);
else

    error(['Unknown platform: ' c]);
    
end
  
if (nargout > 0) | (~isempty(MEX_RETURN_RESULT_IN_ANS)  & ...
        MEX_RETURN_RESULT_IN_ANS)
    errorCode = errCode;
elseif (errCode ~= 0)
    error('Unable to complete successfully');
end


%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function result = read_response_file(filename)
%
% Read a response file (a filename that starts with '@')
% and return a cell of strings, one per entry in the response file.
% Use Perl to ensure processing of arguments is the same as mex.bat
%

result = {};

cmd = ['"' matlabroot '\sys\perl\win32\bin\perl" -e "' ...
    'require ''' matlabroot '\\sys\\perl\\win32\\lib\\shellwords.pl'';' ...
    'open(FILE, ''' filename ''') || die ''Could not open ' filename ''';' ...
    'while (<FILE>) {$line .= $_;} ' ...
    '$line =~ s/\\/\\\\/g;' ...
    '@ARGS = &shellwords($line); ' ...
    '$\" = \"\n\";' ...
    'print \"@ARGS\";'];

[s, r] = dos(cmd);

if s == 0
    cr = sprintf('\n');
    while ~isempty(r)
        [result{end+1}, r] = strtok(r, cr);
    end
end

function [mexname, setup] = get_mex_opts(varargin)
%
% GET_MEX_OPTS gets the options from the command line.
%
% name:
% It gets the name of the destination MEX-file.  This has two
% purposes: 
%   1) All platforms need to clear the MEX-file from memory before
%      attempting the build, to avoid problems rebuilding shared
%      libraries that the OS considers "in use".
%   2) Windows MATLAB deletes the MEX-file before the build occurs.
%      It then checks to see whether the MEX-file was created so as
%      to establish error status.
%   This function returns the minimum necessary information.  Further
%   processing is done on the MEX-file name by clear_mex_file to 
%   successfully clear it.
%
% setup:
% It also returns whether or not '-setup' was passed.
%

mexname = '';
outdir = '';
setup = 0;

% First, check for and expand response files into varargin.
v = {};
for count=1:nargin
    arg = varargin{count};
    if arg(1) == '@'
        new_args = read_response_file(arg(2:end));
        v(end+1:end+length(new_args)) = new_args;          
    else
        v{end+1} = arg;
    end
end

varargin = v;

count = 1;
while (count <= nargin)
    arg = varargin{count};
    if isempty(mexname) & arg(1) ~= '-' & ~any(arg=='=') & any(arg=='.')
        %
        % Source file: MEX-file will be built in current directory
        % Only the first source file matters
        %
        mexname = arg;
        fileseps = find(mexname == filesep);
        if any(fileseps)
            mexname = mexname(fileseps(end)+1:end);
        end
        mexname = strtok(mexname, '.');
    elseif strcmp(arg, '-f')
        count = count + 1;
    elseif strcmp(arg, '-output')
        count = count + 1;
        if count > length(varargin)
            error('-output switch must be followed by a filename');
        end
        mexname = varargin{count};
    elseif strcmp(arg, '-outdir')
        count = count + 1;
        if count > length(varargin)
            error('-outdir switch must be followed by a directory name');
        end
        outdir = varargin{count};
    elseif strcmp(arg, '-setup')
        setup = 1;
        break;
    end
    count = count + 1;
end

mexname = fullfile(outdir, mexname);

function clear_mex_file(basename)
%
% CLEAR_MEX_FILE Clear a MEX-file from memory.  This is a tricky
%   business and should be avoided if possible.  It takes a relative
%   or absolute filename as the MEX-file name, and the list of loaded
%   MEX-file names.
%
%   If CLEAR_MEX_FILE is unable to clear the MEX-file, it will error.
%   This can happen if the MEX-file is locked.

seps = find(basename == filesep);
if isempty(seps)
    % basename is in the current directory
    fullname = fullfile(cd,basename);
else
    % -output was used to determine the location, as well as the
    % name, of the mex file.  The easiest way to find the dir.
    % it's in is to cd to it.
    savedir = cd;
    destdir = basename(1:seps(end));
    cd(destdir);
    fullname = fullfile(cd,basename);
    cd(savedir);
end

if ~isempty(findstr(fullname, 'private'))
    % Things in private directories are represented by the full
    % path
    mexname = fullname;
else
    modifiers = find(fullname == '@');
    if any(modifiers)
        % Methods have the class directory prepended
        mexname = fullname((modifiers(end)+1):end);
        % Methods are always displayed with UNIX file
        % separators
        mexname(mexname==filesep) = '/';
    else
        % Otherwise, we just use the base name
        mexname = basename;
    end
end

clear_mex(mexname);
% Make sure that the MEX-file is cleared
[ms, mexs] = inmem;
if ~isempty(strmatch(mexname, mexs, 'exact'))
    error('Your MEX-file is locked and must be unlocked before recompiling.');
end

function clear_mex(varargin)
% This will clear a MEX-file successfully, because it has no internal
% variables.  varargin is a builtin function and is therefore not a
% valid MEX-file name.

clear(varargin{:});

function [cmd, rspfile] = make_rsp_file(cmdtool, cmdargs)
rspfile = [tempname '.rsp'];
[Frsp, errmsg] = fopen(rspfile, 'wt');
if Frsp == -1
    error('Cannot open file "%s" for writing: %s.', rspfile, errmsg)
end
try
    count = fprintf(Frsp, '%s', cmdargs);
    if count < length(cmdargs)
        errmsg = ferror(Frsp);
        error('Cannot write to file "%s": %s.', rspfile, errmsg);
    end
    fclose(Frsp);
catch
    fclose(Frsp);
    delete(rspfile);
    error(lasterr);
end

    cmd = [cmdtool ' @"' rspfile '"'];

