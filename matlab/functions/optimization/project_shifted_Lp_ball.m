%% project onto the shifted Lp ball
%{
    lp-ball:  || x - shift_vec ||_p^p = n / (2^p)
    shift_x = x-shift_vec; 
    then 
    shift_xp = shift_x / ( ||shift_x||_p / ( n^(1/p)/2 )  ); 
    xp = shift_xp + shift_vec; 
%}
function [xp] = project_shifted_Lp_ball(x,shift_vec, p)
shift_x = x-shift_vec;
normp_shift = norm(shift_x, p);
n = numel(x);

if normp_shift^p~=n/(2^p)
    xp = shift_x / (normp_shift / (n^(1/p)/2) ) + shift_vec;
else
    xp = x;
end

return;
