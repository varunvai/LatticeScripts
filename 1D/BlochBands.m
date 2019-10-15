clear all

% Lattice depth (photon recoil energies)
V0 = 15;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
q = -1:0.01:1; 
Nq = length(q);
Nbase = 21;
En = zeros(Nbase,Nq);

for j = 1:Nq
    H = HLattice(q(j),V0,Nbase);
    En(:,j) = eig(H);
end

