function filesError = compileDir(dirInput);
% compiles mex files in specified directory
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin < 1
    dirInput = pwd;
end

footpath = cd;
addpath(genpath([footpath '/util']));

files = dir2(dirInput);
filesError = compileFiles({files.filepath});
