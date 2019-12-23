function t = fwStepOpt(a, b)
% Compute the optimal step size based on the shape of the quadratic function.
%
% Input
%   a       -  second-order coefficient
%   b       -  first-order coefficient
%           
% Output    
%   t       -  optimal step
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-16-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% linear not parabola
if abs(a) < eps
    if b > 0
        t = 1;
    else
        t = 0;
    end
    return;
end

t = -b / a / 2;
if t <= 0
    if a > 0
        t = 1;
    else
        t = 0;
    end

elseif t <= .5
    if a > 0
        t = 1;
    end
    
elseif t <= 1
    if a > 0
        t = 0;
    end
else
    if a > 0
        t = 0;
    else
        t = 1;
    end
end
