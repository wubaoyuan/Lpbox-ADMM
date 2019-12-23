function [user,host,userAthost] = getUserHost;
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

host='';
user='';

if ispc
    try
        [s,user] = system('ECHO %USERNAME%');
        user(end) = [];
        user = lower(user);
    catch
        user = getUserHostFromDir;
    end
    %     [s,host] = system('ECHO %USERDOMAIN%');
    %     host = lower(host);
    [s,host] = system('hostname');
    host(end) = [];

    %     cygwin = 'C:\tim\applications\cygwin_1.5.10-3\installation_files\bin\';
    %     if ~exist(cygwin)
    %         user = getUserHostFromDir;
    %     else
    %         [s,user] = system([cygwin 'whoami']);
    %         user(end) = [];
    %     end
    %     [s,host] = system('hostname');
    %     host(end) = [];
else
    [s,user] = system('whoami');
    user(end) = [];
    [s,host] = system('hostname');
    host(end) = [];
end


if isempty(host) 
    userAthost = user;
else
    userAthost = [user '@' host];
end

u1 = regexp(userAthost, 'tcour@v0[1-9]-64.vision.grasp.upenn.edu');
u2 = regexp(userAthost, 'tcour@v0[1-9]-64');
u3 = regexp(userAthost, 'tcour@v0[1-9]-32.vision.grasp.upenn.edu');
u4 = regexp(userAthost, 'tcour@v0[1-9]-32');
if length(u1)||length(u2)||length(u3)||length(u4)|| strcmp(userAthost,'tcour@super.vision.grasp.upenn.edu')
    userAthost = 'tcour@vision.grasp.upenn.edu';
end

%TODO:getenv('HOSTNAME')

function user = getUserHostFromDir;
if exist('D:/research/projets Matlab/projet graphcut/graphcut','dir')==7
    user = 'tim labo';
elseif exist('C:\gsong\Source\memory_graph\corel_img\selected\www\','dir')==7
    user = 'gsong laptop';
else
    error('Error : user not recognized. Add your username, change your home, image, and results directories');
end


