function gphs = gphEgNewAsg(asg, gph0s)
% Keep only edge with correspondence.
%
% Input
%   asg     -  assignment
%   gph0s   -  original graph, 1 x 2 (cell)
%
% Output
%   gphs    -  new graph, 1 x 2 (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-20-2012

gph10 = gph0s{1};
gph20 = gph0s{2};

% edge
Eg10 = gph10.Eg;
Eg20 = gph20.Eg;

% dimension
X = asg.X;
[n1, n2] = size(X);
m1 = size(Eg10, 2);
m2 = size(Eg20, 2);

% index of matched points
idx = find(X(:));
[i1s, i2s] = ind2sub([n1, n2], idx);

% remove edge
Eg1 = removeEg(Eg10, i1s);
Eg2 = removeEg(Eg20, i2s);

% store
gph1 = gph10;
gph1.Eg = Eg1;
gph2 = gph20;
gph2.Eg = Eg2;

gphs = {gph1, gph2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Eg = removeEg(Eg0, is)
% remove edge

% dimension
m = size(Eg0, 2);

vis = zeros(1, m);
for c = 1 : m
    if any(Eg0(1, c) == is) && any(Eg0(2, c) == is)
        vis(c) = 1;
    end
end

Eg = Eg0(:, vis == 1);
