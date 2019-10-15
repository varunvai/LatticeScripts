function [U] = FreeEvolve(k,psi_0,t)

H = HLattice(k,0,length(psi_0));
U = expm(-i*2*pi*H.*t);

end

