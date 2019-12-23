function [x] = prox_box_r01(a,b,c)
% solve the following OP:
% min_{x} 0.5 ||x-a||^2 + <|x|,b> + <|x-1|,c>, s.t. 0 <= x <= 1
% we assume that both b and c are non-negative

n = length(a);
x = zeros(n,1);
for i=1:n,
    x1 = 0;
    x2 = 1;
    x3 = max(0,min(1,a(i)+c(i)-b(i)));
    f1 = 0.5*a(i)*a(i)+c(i);
    f2 = 0.5*(a(i)-1)^2 + b(i);
    f3 = 0.5*(x3-a(i))^2 + b(i)*x3 + (1-x3)*c(i);
    if(f1<=f2 && f1<=f3)
        x(i) = x1;
    elseif(f2<=f1 && f2<=f3)
        x(i) = x2;
    elseif(f3<=f1 && f3<=f2)
        x(i) = x3;
    end
end


