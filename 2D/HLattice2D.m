function [H,M,N] = HLattice2D(qx,qy,V0,Nbase)
%%   HLattice2D Returns the hamiltonian for the RbYb 2D lattice
%   This lattice has 3 beams with k-vectors along x, -x and y.
%   This code assumes equal intensities for all beams.
%   Will be fixed to have arbitrary intensities later
%
%   Inputs:
%   qx,qy = Quasimomentum in recoil momenta
%   V0 = Lattice "depth" in recoil energies
%   Nbase = # of k vectors in each direction
%
%   Outputs:
%   H = The Hamiltonian
%   M = Plane wave index along x
%   N = Plane wave index along y

%% Set up indexing for basis vectors
M = -Nbase:Nbase;
N = M;
for j = 1:Nbase
    M = [(-Nbase:Nbase)-j M (-Nbase:Nbase)+j];
    N = [(-Nbase:Nbase)+j N (-Nbase:Nbase)-j];
end

%% Define lattice hamiltonian
H = zeros((2*Nbase+1)^2,(2*Nbase+1)^2);

% Determine plane-wave couplings for the RbYb 2D lattice
for j = 1:length(M)
    for l = 1:length(M)
        if (M(l) == M(j)+2) & (N(l) == N(j))
            H(j,l) = -V0;
        elseif (M(l) == M(j)-1) & (N(l) == N(j)-1)
            H(j,l) = -V0;
        elseif (M(l) == M(j)+1) & (N(l) == N(j)-1)
            H(j,l) = -V0;
        end
    end
end

% Make sure the Hamiltonian is Hermitian and set diagonal elements
H = H + H';
H = H + diag((qx+M).^2 + (qy+N).^2-3*V0);

end

