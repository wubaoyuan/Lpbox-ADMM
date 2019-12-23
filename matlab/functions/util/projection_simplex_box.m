function [x] = projection_simplex_box(a,b,c,k)
% This program solves the following OP:
% min_{x} 0.5 |x-a|_2^2
% s.t. 0 <= x <= b, x'c >= k
% Warning:
% We assume b>=0, c>=0 and b'c>=k

debug = 0;
if(debug && k>sum(b.*c)),
    error('illegal k!');
end

% Any entry below eps1 will be set to 0 automatically
eps1 = eps;
n = length(a);
low = inf;
up = -inf;

for i=1:n,
    ac = a(i) / c(i);
    bc = b(i) / c(i);
    if(c(i) > eps1)
        low = min(low,-abs(ac-bc));
        up = max(up, ac);
    end
end

low = low - 1;
up = up + 1;

if(debug)
    check_lowerbound(low,a,b,c);
    check_upperbound(up,a,b,c);
end

fobj_handle = @ (lambda) c'* max(0,min(b,a-lambda*c)) - k;
lambda = bisection(fobj_handle, low, up);
x = max(0,min(b,a-lambda*c));

function [flag]= check_lowerbound(lambda,a,b,c)
% min(b,a-lambda c) = b, for all ci ~=0
n = length(a);
flag = 0;
eps1 = eps;
for i=1:n,
    if (abs(c(i))>eps1 && a(i)-lambda*c(i) < b(i) ),
        flag = 1;
        lambda
        a(i)
        b(i)
        c(i)
        a(i)-lambda*c(i) - b(i)
        error('wrong!');
        break;
    end
end


function [flag]= check_upperbound(lambda,a,b,c)
% min(b,a-lambda c) <= 0
% a-lambda c <= 0
n = length(a);
flag = 0;
eps1 = eps;
for i=1:n,
    if( c(i) > eps1  && a(i)-lambda*c(i)>0 ),
        a(i)
        a(c)
        a(i)-lambda*c(i)
        flag = 1;
        error('wrong!');
        break;
    end
end


function [ret,iter] = bisection( f, a, b)

eps1 = 1e-12;
a1 = a;
b1 = b;
fa = f(a);
fb = f(b);
if ( fa == 0 )
    ret = a;
    return;
elseif ( fb == 0 )
    ret = b;
    return;
elseif ( fa * fb > 0 )
       ret = (a+b)/2;
%     error( 'f(a) and f(b) do not have opposite signs' );
return;
end

max_iter = 300;
for iter = 1:max_iter,
    c = (a + b)/2; 
    fc = f(c);
    if (abs(fc)<=eps1), ret = c; break; end
    fa = f(a);
    if (fc*fa < 0 )
        b = c;
    else
        a = c;
    end
end

if(max_iter==iter),
% fprintf('maximum iteration! iter:%d, f:%.5e, a:%.5e, b:%.5e\n',iter,abs(fc),a1,b1);
ret = c;
end

