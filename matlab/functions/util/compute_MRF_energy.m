function [c] = compute_MRF_energy(x,A,b)
c = x'*A*x+b'*x;
return;