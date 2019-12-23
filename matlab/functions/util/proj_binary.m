function [x] = proj_binary(a)
% solve the following QP:
% min_{x} 0.5 ||x-a||^2, s,t, x\in{0,1}

n = length(a);
x = zeros(n,1);
for i=1:n,
    if(a(i)>0.5)
        x(i)=1;
    else
        x(i)=0;
    end
end


