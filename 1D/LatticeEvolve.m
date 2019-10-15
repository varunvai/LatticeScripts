function [U] = LatticeEvolve(HLatt,psi_0,t)

[V,En] = eig(HLatt);
U = V*expm(-i*2*pi*En.*t)*V';

end
