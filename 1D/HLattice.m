function [H] = HLattice(k,V0,Nbase)
%HLattice Returns the hamiltonian for a lattice
%   Inputs:
%   k = Quasimomentum in recoil momenta
%   V0 = Lattice depth in recoil energies
%   Nbase = # of Brillouin zones to sample
%
%   Outputs:
%   H = The Hamiltonian

H = zeros(Nbase,Nbase);

for i=1:Nbase
    H(i,i)=(k+2*(i-(Nbase+1)/2)).^2;
end

H = H+diag(-0.25*V0*ones(Nbase-1,1),1)+diag(-0.25*V0*ones(Nbase-1,1),-1);

end

