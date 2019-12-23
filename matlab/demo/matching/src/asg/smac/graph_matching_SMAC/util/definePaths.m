function [paths,user,host,userAthost] = definePaths(user,host);
%{
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

global data;
[data.paths,data.user,data.host] = definePaths;
%}

if nargin == 0
    [user,host,userAthost] = getUserHost;
end
if strcmp(user,'timothee') && ismac
    host='MacBookPro.local';
    userAthost='timothee@MacBookPro.local';    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% static directories to change to your own directories %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0 % means that paths are relative to local machine
    paths.home = pwd;
end

switch userAthost
    case 'timothee@MacBookPro.local'
        paths.java='/usr/bin/java';
        
        paths.dir_parser_stanford='/Users/timothee/research/programs/autre_code/text_processing/stanford-parser-2007-08-19';
        paths.dir_parser='/Users/timothee/research/programs/autre_code/text_processing/eleni';
        paths.dirInput_movies='/Users/timothee/research/data/video/movies';
%         paths.dirInput_movies='/Volumes/MyBook/research/data/movies/';%external
%         paths.dir_screenplays='/Users/timothee/research/data/text/movie scripts/data/files';


        
        %data
        paths.repertoireImages = '/Users/timothee/research/data/images';
        paths.fichierImageInitiale = '/Users/timothee/research/data/images/108005.jpg';
        paths.videoDirectory = '/Users/timothee/research/data/videos';
        paths.repertoireResults = '/Users/timothee/research/results';
        paths.repertoireBenchmarks = '/Users/timothee/research/results/benchmarks';
        paths.savedWorkspace = '/Users/timothee/research/results/savedWorkspace';

        paths.nistTraining = '/Users/timothee/research/data/images/NIST/nistTraining.mat';
        paths.nistTesting = '/Users/timothee/research/data/images/NIST/nistTesting.mat';

%         paths.dir_mplayer='/Users/timothee/research/programs/mplayer';
        paths.dir_mplayer='/Users/timothee/programs/mplayer';
        
        %other
        paths.oldFilesDirectory = '/Users/timothee/research/programs/old graphcuts/';
%         paths.temp = fullfile(paths.savedWorkspace,'temp');
        paths.dir_temp='/Users/timothee/temp/trash';
        
        %data
        paths.dir_hand_labels='/Users/timothee/research/data/hand labels';
        paths.dir_segmented_database='/Users/timothee/research/results/_segmented databases';
        paths.dir_horse_image='/Users/timothee/research/data/images/object categories/images with ground truth/horses/image';
        paths.dir_horse_mask='/Users/timothee/research/data/images/object categories/images with ground truth/horses/mask';
        paths.dir_horse_classes='/Users/timothee/research/results/_segmented databases/horses';        
    case 'timothee@M400'
        %data
        paths.repertoireImages = 'C:\tim\research\data\images';
        paths.fichierImageInitiale = 'C:\tim\research\data\images\108005.jpg';
        paths.videoDirectory = 'C:\tim\research\data\videos';
        paths.repertoireResults = 'C:\tim\research\results';
        paths.repertoireBenchmarks = 'C:\tim\research\results\benchmarks';
        paths.savedWorkspace = 'C:\tim\research\results\savedWorkspace';

        paths.nistTraining = 'C:\tim\research\data\images\NIST\nistTraining.mat';
        paths.nistTesting = 'C:\tim\research\data\images\NIST\nistTesting.mat';

        paths.cygwin = 'C:\programs\cygwin\bin';
        paths.sshFile = 'C:\programs\cygwin\bin\ssh';
        paths.scpFile = 'C:\programs\cygwin\bin\scp -c blowfish -C';
        paths.msPaint = 'C:\WINDOWS\System32\mspaint.exe';
        paths.jpeg2ps = 'C:\programs\bin\jpeg2ps-1.9\jpeg2ps';
        paths.irfanview = '"C:\Program Files\IrfanView\i_view32.exe"';
        paths.imageMagick = '"C:\Program Files\ImageMagick-6.2.7-Q16\convert.exe"';
        paths.graphicsMagick = '"C:\Program Files\GraphicsMagick-1.1.7-Q8\gm.exe"';

        %other
        paths.oldFilesDirectory = 'C:\tim\research\programs\old graphcuts\';
        paths.temp = fullfile(paths.savedWorkspace,'temp');
        
        %data
        paths.dir_hand_labels='C:\tim\research\data\hand labels';
        paths.dir_segmented_database='C:\tim\research\results\_segmented databases';
        paths.dir_horse_image='C:\tim\research\data\images\object categories\images with ground truth\horses\image';
        paths.dir_horse_mask='C:\tim\research\data\images\object categories\images with ground truth\horses\mask';
        paths.dir_horse_classes='C:\tim\research\results\_segmented databases\horses';
    case {'timothee@alef.seas.upenn.edu','timothee@bet.seas.upenn.edu','timothee@gimel.seas.upenn.edu','timothee@dalet.seas.upenn.edu',}
        paths.java='/usr/java/jdk1.6.0/bin/java';
        paths.dir_parser_stanford='~/programs/toolbox/text_processing/stanford-parser-2007-08-19';
        paths.dir_parser='~/programs/toolbox/text_processing/eleni';

        paths.dirInput_movies='/mnt/zain/exports/users/timothee/movies/lost';
        paths.dir_screenplays='/mnt/zain/exports/users/timothee/movies/screenplays/data/files';%VOIR
        
    case 'tcour@vision.grasp.upenn.edu'
        paths.graphcut = '/home/tcour/programs/graphcut';
        paths.repertoireResults = '/home/tcour/results/';
        paths.savedWorkspace = '/home/tcour/savedWorkspace/';
        paths.nistTraining = '/home/tcour/data/NIST/nistTraining.mat';
        paths.nistTesting = '/home/tcour/data/NIST/nistTesting.mat';
        paths.repertoireImages = 'data/images/';
        paths.fichierImageInitiale = 'data/images/baby.jpg';
        paths.oldFilesDirectory = '/home/tcour/programs/old graphcuts/';
        
        %data
        paths.dir_segmented_database='/data/insecure/tim/results/segmentations/';
        paths.dir_horse_mask='/home/tcour/data/images with ground truth/horses/mask';
        paths.dir_horse_classes='/data/insecure/tim/results/segmentations/horses_new';

        paths.dirInput_movies='/data/insecure/tim/videscribe/movies';
        paths.dir_screenplays='/data/insecure/tim/videscribe/screenplays/data/files';
        
        
    case 'timothee@timothee-dell'
        %data
        paths.repertoireImages = 'C:\tim\research\data\images';
        paths.fichierImageInitiale = 'C:\tim\research\data\images\108005.jpg';
        paths.videoDirectory = 'C:\tim\research\data\videos';
        paths.repertoireResults = 'C:\tim\research\results';
        paths.repertoireBenchmarks = 'C:\tim\research\results\benchmarks';
        paths.savedWorkspace = 'C:\tim\research\results\savedWorkspace';
        
        paths.cygwin = 'C:\programs\cygwin\bin';
        paths.sshFile = 'C:\programs\cygwin\bin\ssh';
        paths.scpFile = 'C:\programs\cygwin\bin\scp -c blowfish -C';
        paths.msPaint = 'C:\WINDOWS\System32\mspaint.exe';
        paths.jpeg2ps = 'C:\programs\bin\jpeg2ps-1.9\jpeg2ps';

        %other
        paths.oldFilesDirectory = 'C:\tim\research\programs\old graphcuts\';
        paths.temp = fullfile(paths.savedWorkspace,'temp');
        
    case 'administrator@FUJITSU'
        %data
        paths.repertoireImages = 'C:\tim\research\data\images\chinese database\images\mouth';
        paths.fichierImageInitiale = 'C:\tim\research\data\images\chinese database\images\mouth\1.jpg';
        paths.videoDirectory = 'C:\tim\research\data\videos';
        paths.repertoireResults = 'C:\tim\research\results';
        paths.repertoireBenchmarks = 'C:\tim\research\results\benchmarks';

        %workspace
        paths.savedWorkspace = 'C:\tim\research\saved worspace files';
        paths.temp = 'C:\tim\research\saved worspace files\temp';

        %databases
        paths.nistTraining = 'C:\tim\research\data\images\nist\nistTraining.mat';
        paths.nistTesting = 'C:\tim\research\data\images\nist\nistTesting.mat';
        %paths.cmuFaces = 'C:\tim\programmation\data\images\package faces vs nonfaces\data\faces.pat';
        %paths.cmuNonFaces = 'C:\tim\programmation\data\images\package faces vs nonfaces\data\nonfaces.pat';

        %programs
%         paths.irfanview = '"c:\program files\irfanview2\i_view32"';
        paths.irfanview = '"C:\Program Files\IrfanView\i_view32.exe"';
        paths.msPaint = 'C:\WINDOWS\System32\mspaint.exe';
        paths.jpeg2ps = 'C:\programs\bin\jpeg2ps-1.9\jpeg2ps';
        %paths.eps2jpg = 'C:\Program Files\autre\FigEpsPdf2004-05-23\eps2jpg';
        paths.notpad = 'C:\windows\notpad.exe';

        paths.cygwin = 'C:\programs\cygwin\bin';
        paths.sshFile = 'C:\programs\cygwin\bin\ssh';
        paths.scpFile = 'C:\programs\cygwin\bin\scp -c blowfish -C';

        %other
        paths.oldFilesDirectory = 'C:\tim\research\programs\old graphcuts\';

    case 'timothee@ba'
        %data
        paths.repertoireImages = 'C:/tim/programmation/data/images/ma base/';
        paths.fichierImageInitiale = 'C:/tim/programmation/data/images/ma base/14.jpg';
        paths.videoDirectory = 'C:/tim/programmation/data/videos/';
        paths.repertoireResults = 'C:/tim/programmation/projetGraphcut/results/';
        paths.repertoireBenchmarks = 'C:\tim\programmation\projetGraphcut\results\benchmarks\';

        %workspace
        paths.savedWorkspace = 'C:\tim\programmation\projetGraphcut\saved worspace files\';
        paths.temp = 'C:\tim\programmation\projetGraphcut\saved worspace files\temp\';

        %databases
        paths.nistTraining = 'C:\tim\programmation\data\images\NIST\nistTraining.mat';
        paths.nistTesting = 'C:\tim\programmation\data\images\NIST\nistTesting.mat';
        paths.cmuFaces = 'C:\tim\programmation\data\images\package faces vs nonfaces\data\faces.pat';
        paths.cmuNonFaces = 'C:\tim\programmation\data\images\package faces vs nonfaces\data\nonfaces.pat';

        %programs
        paths.irfanview = '"c:\program files\irfanview2\i_view32"';
        paths.msPaint = 'C:\WINDOWS\System32\mspaint.exe';
        paths.jpeg2ps = 'C:\Program Files\autre\jpeg2ps-1.9\jpeg2ps';
        %paths.eps2jpg = 'C:\Program Files\autre\FigEpsPdf2004-05-23\eps2jpg';
        paths.notpad = 'C:\windows\notpad.exe';

        paths.cygwin = 'C:\tim\applications\cygwin_1.5.10-3\installation_files\bin\';
        paths.sshFile = 'C:\tim\applications\cygwin_1.5.10-3\installation_files\bin\ssh';
        paths.scpFile = 'C:\tim\applications\cygwin_1.5.10-3\installation_files\bin\scp -c blowfish -C';

        %other
        paths.oldFilesDirectory = 'C:\tim\programmation\projetGraphcut\old graphcuts\';

    case 'tim labo@timothee' %????????????
        paths.repertoireImages = 'D:/research/images/ma base';
        paths.fichierImageInitiale = 'D:/research/images/ma base/14.jpg';
        paths.repertoireResults = 'D:/research/projets Matlab/projet graphcut/results';
        paths.jpeg2ps = '?';
        paths.nistTraining = 'D:\research\projets Matlab\projet graphcut\saved workspace\NIST\nistTraining.mat';
        paths.nistTesting = 'D:\research\projets Matlab\projet graphcut\saved workspace\NIST\nistTesting.mat';
        paths.oldFilesDirectory = 'D:\research\projets Matlab\projet graphcut\old graphcuts\';
    case 'matlab@lvn403-17.grasp'%.upenn.edu'

        homeDir = '/home/matlab/tim';
        paths.repertoireImages = '/home/matlab/tim/data/ma base';
        paths.fichierImageInitiale = '/home/matlab/tim/data/ma base/14.jpg';
        paths.repertoireResults = '/home/matlab/tim/results';
        %paths.jpeg2ps = '?';
        paths.cmuFaces = '/home/matlab/tim/data/faces.pat';
        paths.cmuNonFaces = '/home/matlab/tim/data/nonfaces.pat';
        paths.nistTraining = '/home/matlab/tim/data/nistTraining.mat';
        paths.nistTesting = '/home/matlab/tim/data/nistTesting.mat';
        paths.oldFilesDirectory = '/home/matlab/tim/programs/old graphcuts';
    case 'jianbo home' %?? '/home/jshi/digit/@??'


        userName = 'jianbo home';
        paths.repertoireImages = '/home/jshi/digit/';
        paths.fichierImageInitiale = '/home/jshi/digit/I_base_dives.jpg';
        paths.repertoireResults = '/home/jshi/digit/';
        paths.jpeg2ps = '?';
    case 'jianbo labo' %? '/home/matlab/jshi/digit/'
        paths.repertoireImages = '/home/matlab/jshi/digit/';
        paths.fichierImageInitiale = '/home/matlab/jshi/digit/I_base_dives.jpg';
        paths.repertoireResults = '/home/matlab/jshi/digit/';
        paths.jpeg2ps = '?';

    case 'gsong laptop' %? 'C:\gsong\Source\memory_graph\corel_img\selected\www\'
        paths.repertoireImages = 'C:\gsong\Source\memory_graph\corel_img\selected\www\';
        paths.fichierImageInitiale = 'C:\gsong\Source\memory_graph\corel_img\selected\www\50mb_rectangle_business_card_cd-rom_silk-screened.jpg';
        paths.repertoireResults = 'C:\gsong\Source\memory_graph\corel_img\selected\www\';
        paths.cmuFaces = 'C:\gsong\Source\gang\data\data\faces.pat';
        paths.cmuNonFaces = 'C:\gsong\Source\gang\data\data\nonfaces.pat';
        paths.jpeg2ps = '?';
        
    case 't-ticour@t-ticour'
        paths.graphcut = 'C:\timothee\research\code\graphcut';
        paths.repertoireResults = 'C:\timothee\research\results\';
        paths.savedWorkspace = 'C:\timothee\research\savedWorkspace\';
        paths.temp = 'C:\timothee\research\temp\';
        
        paths.nistTraining = 'C:\timothee\research\data\NIST\nistTraining.mat';
        paths.nistTesting = 'C:\timothee\research\data\NIST\nistTesting.mat';
        paths.repertoireImages = 'C:\timothee\research\data\ma base\';
        paths.fichierImageInitiale = 'C:\timothee\research\data\ma base\14.jpg';
        paths.oldFilesDirectory = 'C:\timothee\research\code\old graphcuts';

        paths.irfanview = '"c:\program files\irfanview\i_view32"';
        
    otherwise
        disp('You can change your home, image, and results directories if you want ; see startup/definePaths');
%         paths.repertoireImages = 'data/images/';
%         paths.fichierImageInitiale = 'data/images/baby.jpg';
%         paths.repertoireResults = '';
end

paths.includeDir = fullfile(pwd,'include');
