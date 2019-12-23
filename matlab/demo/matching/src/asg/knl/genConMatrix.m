function C = genConMatrix(seg1, seg2)
% Generate the confusion matrix of two segmentations.
%
% Input
%   seg1    -  1st segmentation
%   seg2    -  2nd segmentation
%
% Output
%   C       -  confusion matrix, k1 x k2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-28-2012

% segmentation
s1 = seg1.s; G1 = seg1.G; [k1, m1] = size(G1); l1 = G2L(G1);
s2 = seg2.s; G2 = seg2.G; [k2, m2] = size(G2); l2 = G2L(G2);

% overlap of each pair of segments
O = zeros(m1, m2);
for i1 = 1 : m1
    for i2 = 1 : m2
        a = max(s1(i1), s2(i2));
        b = min(s1(i1 + 1), s2(i2 + 1));
        if a < b
            O(i1, i2) = b - a;
        end
    end
end

% confusion matrix
C = zeros(k1, k2);
for c1 = 1 : k1
    vis1 = l1 == c1;
    for c2 = 1 : k2
        vis2 = l2 == c2;

        C(c1, c2) = sum(sum(O(vis1, vis2)));
    end
end

% I1 = seg2IH(seg1);
% I2 = seg2IH(seg2);
% C = I1 * I2';
% equal('C', C, C2);
