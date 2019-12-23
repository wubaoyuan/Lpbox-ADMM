function [x,y] = co2xy(z,m)
y = floor(z/m);
t = mod(z,m);
if(t>0)
    y = y+1;x = t;
else
    x = m;
end



