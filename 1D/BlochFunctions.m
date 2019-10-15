clear all

% Lattice depth (photon recoil energies)
V0 = 5;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
q = 1; band = 1;
x = -2*pi:0.1:2*pi;
Nx = length(x);
Nbase = 21;
n = -(Nbase-1)/2:(Nbase-1)/2;
phi = x;

for j = 1:Nx
    H = HLattice(q,V0,Nbase);
    [state, En] = eig(H);
    phase = exp(i*2*n*x(j));
    phi(j) = phase*state(:,band);
            
end
