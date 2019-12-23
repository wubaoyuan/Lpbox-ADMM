function H = mconv(alg, F, g, l)
% Matrix convolution.
%
% Input
%   alg     -  algorithm type, 1 | 2 | 3 | 4
%              1: normal
%              2: truncated
%              3: truncated
%   F       -  template, d x m
%   g       -  impulse, n x 1
%   l       -  parameter, used only when alg == 2 or 3
% 
% Output
%   H       -  result, d x length
%              length = (n + m - 1), if alg == 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-10-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
[d, m] = size(F);
n = size(g, 1);

% flip
if alg == 2 || alg == 3 || alg == 4
    g = g(end : - 1 : 1);
end

% 1-D convolution in each dimension
H = convFast(F, g');
% H2 = convFast2(F, g');
% equal('H', H, H2);

% normal
if alg == 1

% truncated 
elseif alg == 2
    tail = min(n + m - 1, n + l - 1);
    H = H(:, n : tail);

% truncated
elseif alg == 3
    head = max(1, n - l + 1);
    tail = min(n + m - 1, n + l - 1);
    H = H(:, head : tail);
    
% truncated 
elseif alg == 4
    n = round((n + 1) / 2);
    tail = min(n + m - 1, n + l - 1);
    H = H(:, n : tail);

else
    error('unknown algorithm: %d', alg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = convFast(F, g)
% Fast convolution of two vectors.

% dimension
[d, m] = size(F);
n = size(g, 2);

Ly = m + n - 1;

% Find smallest power of 2 that is > Ly
global Ly2s;
if isempty(Ly2s)
    Ly2 = pow2(nextpow2(Ly));
    Ly2s = [Ly; Ly2];
else
    isFound = 0;
    for i = 1 : size(Ly2s, 2)
        if Ly2s(1, i) == Ly
            Ly2 = Ly2s(2, i);
            isFound = 1;
            break;
        end
    end
    if isFound == 0
        Ly2 = pow2(nextpow2(Ly));
        Ly2s = [Ly2s, [Ly; Ly2]];
    end
end

% Fast Fourier transform
G = fft(g, Ly2);

H = zeros(d, n + m - 1);
for c = 1 : d
    % Fast Fourier transform
    F2 = fft(F(c, :), Ly2);

    Y = F2 .* G;

    % Inverse fast Fourier transform
    y = real(ifft(Y, Ly2));

    H(c, :) = y(1 : 1 : Ly);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = convFast2(F, g)
% Fast convolution of two vectors.

% dimension
[d, m] = size(F);
n = size(g, 2);

Ly = m + n - 1;

% Find smallest power of 2 that is > Ly
Ly2 = pow2(nextpow2(Ly));

% Fast Fourier transform
G = fft(g, Ly2);

H = zeros(d, n + m - 1);
for c = 1 : d
    % Fast Fourier transform
    F2 = fft(F(c, :), Ly2);

    Y = F2 .* G;

    % Inverse fast Fourier transform
    y = real(ifft(Y, Ly2));

    H(c, :) = y(1 : 1 : Ly);
end
