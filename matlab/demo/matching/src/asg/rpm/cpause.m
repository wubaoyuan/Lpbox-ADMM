function return_value = cpause
% Authors: Haili Chui and Anand Rangarajan
% Date:    04/27/2000
% 
% Contact Information:
%
% Haili Chui:		chui@noodle.med.yale.edu
% Anand Rangarajan:	anand@noodle.med.yale.edu
% 
% Terms:	  
% 
% The source code (M-files) are provided under the
% terms of the GNU General Public License with an explicit
% clause permitting the execution of the M-files from within
% a MATLAB environment. See the LICENSE file for details.
%
% cpress_key.m
% ------------------------------------------------------------------- 
% Pause and Press a key. 
%
% ------------------------------------------------------------------- 

key_input = input('press any key when ready: ', 's');

return_value = 0;
if strcmp(key_input, 'q')
    return_value = 1;
end
