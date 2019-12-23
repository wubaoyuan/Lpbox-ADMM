function [CD, CD1, CD2] = genConMatrixD1(seg1, seg2)
% Generate the confusion matrix of two segmentations.
% The confusion is defined by the length of frames that are overlapped by the same group of segments.
%
% Input
%   seg1    -  1st segmentation
%     s     -  starting position, 1 x (m1 + 1)
%     G     -  indicator matrix, k1 x m1
%   seg2    -  2nd segmentation
%     s     -  starting position, 1 x (m2 + 1)
%     G     -  indicator matrix, k2 x m2
%
% Output
%   CD      -  confusion matrix (discrete), k1 x k2, #seg1 & #seg2 / #seg1 | #seg2 >= .5
%   CD1     -  confusion matrix (discrete), k1 x k2, #seg1 & #seg2 / #seg1 >= .5
%   CD2     -  confusion matrix (discrete), k1 x k2, #seg1 & #seg2 / #seg2 >= .5
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% I1 = seg2IH(seg1); seg1 = I2Seg(I1);
I2 = seg2IH(seg2); seg2 = I2Seg(I2);

s1 = seg1.s; G1 = seg1.G; [k1, m1] = size(G1); l1 = G2L(G1); ns1 = diff(s1);
s2 = seg2.s; G2 = seg2.G; [k2, m2] = size(G2); l2 = G2L(G2); ns2 = diff(s2);

% overlap between segments
O0 = zeros(m1, m2);
for i1 = 1 : m1
    for i2 = 1 : m2
        a = max(s1(i1), s2(i2));
        b = min(s1(i1 + 1), s2(i2 + 1));
        if a < b
            O0(i1, i2) = b - a;
        end
    end
end

% #seg1 & #seg2 / #seg1 >= .5
N1 = repmat(ns1', [1 m2]);
O1 = O0 ./ N1 >= .5;

% #seg1 & #seg2 / #seg2 >= .5
N2 = repmat(ns2, [m1, 1]);
O2 = O0 ./ N2 >= .5;

% #seg1 & #seg2 / #seg1 | #seg2 >= .5
N = N1 + N2 - O0;
O = O0 ./ N >= .5;

% confusion matrix
[CD, CD1, CD2] = zeross(k1, k2);
for c1 = 1 : k1
    vis1 = l1 == c1;
    for c2 = 1 : k2
        vis2 = l2 == c2;

        CD1(c1, c2) = sum(sum(O1(vis1, vis2)));
        CD2(c1, c2) = sum(sum(O2(vis1, vis2)));
        CD(c1, c2) = sum(sum(O(vis1, vis2)));
    end
end
