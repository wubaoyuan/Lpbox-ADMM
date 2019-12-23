function enabledberror2( hasDBError, dbstat )
%ENABLEDBERROR2 turns on "dbstop if error" if hasDBError is true.
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.4.4.1 $  $Date: 2004/12/20 16:45:12 $

if (hasDBError && ~isempty(dbstat))
    if (any(strcmp({dbstat.cond},'error')))
        dbstop if error;
    end
    if (any(strcmp({dbstat.cond},'caught error')))
        dbstop if caught error;
    end
end

