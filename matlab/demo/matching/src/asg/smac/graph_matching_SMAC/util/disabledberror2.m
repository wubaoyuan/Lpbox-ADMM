function [hasDBError, dbstat] = disabledberror2(allerror)
%DISABLEDBERROR2 turns off "dbstop if error" and returns previous
%status.
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.4.4.1 $  $Date: 2004/12/20 16:45:11 $

if (nargin == 0)
    allerror = false;
end

% turn off "dbstop if error" while in the try-catch
dbstat = dbstatus;
hasDBError = false;

if ~isempty(dbstat)
  if (any(strcmp({dbstat.cond},'error')))
      hasDBError = true;
      dbclear if error;
  end
  if (allerror && any(strcmp({dbstat.cond},'caught error')))
      hasDBError = true;
      dbclear if caught error;
  end
end
