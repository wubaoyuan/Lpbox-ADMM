function iPart = coTch(path, nPart)
% Access the count number which has been stored in the specified path.
%
% Input
%   path    -  count path
%   nPart   -  #total number
%
% Output
%   iPart   -  part index
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-12-2012

if ~exist(path, 'file')
    iPart = 0;

else
    % load
    [iPart, nPart] = matFld(path, 'iPart', 'nPart');
end

% increase ++
iPart = iPart + 1;
if iPart > nPart
    pr('done');
    return;
end
   
% save
pr('co: %d/%d', iPart, nPart);
save(path, 'iPart', 'nPart');
