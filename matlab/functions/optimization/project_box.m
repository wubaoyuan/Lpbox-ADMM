%% project onto the box constraints
function [xp] = project_box(x)
xp = x;
xp(x>1)=1;
xp(x<0)=0;
return;