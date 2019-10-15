function [H,M,N] = HLattice2Dtest(qx,qy,couplings,elems,Nbase)
%%   HLattice2D Returns the hamiltonian for the RbYb 2D lattice
%   This lattice has 3 beams with k-vectors along x, -x and y.
%   This code assumes equal intensities for all beams.
%   Will be fixed to have arbitrary intensities later
%
%   Inputs:
%   qx,qy = Quasimomentum in recoil momenta
%   couplings = kj-ki vectors
%   elems = strength of coupling from kj-ki
%   Nbase = # of k vectors in each direction
%
%   Outputs:
%   H = The Hamiltonian
%   M = Plane wave index along x
%   N = Plane wave index along y

%% Set up indexing for basis vectors
% M is the momentum (in recoil units) in the x-direction
% M is the momentum (in recoil units) in the y-direction
M = -Nbase:Nbase;
N = M;
for j = 1:Nbase
    M = [(-Nbase:Nbase)-j M (-Nbase:Nbase)+j];
    N = [(-Nbase:Nbase)+j N (-Nbase:Nbase)-j];
end

%% Define lattice hamiltonian
H = zeros((2*Nbase+1)^2,(2*Nbase+1)^2);

% Only set elements that are different by 'couplings', to non-zero
for j = 1:length(M)
    for l = 1:length(M)
        [j l];
        diff = repmat([M(l) N(l)]-[M(j) N(j)],size(couplings,1),1);
        check = find((sum((diff == couplings),2) == 2));
        H(j,l) = -elems(check);
    end
end

% Make sure the Hamiltonian is Hermitian and set diagonal elements
H = H + H';
H = H + diag((qx+M).^2 + (qy+N).^2-3);

end

