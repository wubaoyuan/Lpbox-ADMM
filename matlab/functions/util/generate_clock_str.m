%% generates a string from the current clock
function [c] = generate_clock_str()
a = round(vec(clock));
c = sprintf('%d_%d_%d_%d_%d_%d',a(1),a(2),a(3),a(4),a(5),a(6));
return