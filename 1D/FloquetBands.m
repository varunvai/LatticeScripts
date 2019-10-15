clear all

% Lattice depth (photon recoil energies)
V0 = 5;
% Modulation depth (photon wavelengths) and frequency (recoil energies)
Delta = 0.01;
omega = 4.25;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
q = -1:0.01:1; 
Nq = length(q);
Nbase = 5;
En = zeros(3*Nbase,Nq);

for j = 1:Nq
    H = HFloquet(q(j),V0,Nbase,Delta,omega);
    En(:,j) = eig(H);
end

