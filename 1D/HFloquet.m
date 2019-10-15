function [H] = HFloquet(k,V0,Nbase,Delta,omega)
%HLattice returns the Floquet hamiltonion for a 1D shaken lattice
%   k - quasi-momentum
%   V0 - lattice depth
%   Nbase - number of basis vectors (must be odd)
%   Delta - modulation depth (fraction of light wavelength)
%   omega - modulation frequency (recoil energies)

Veff = besselj(0,2*2*pi*Delta)*V0;
Hbloch = HLattice(k,Veff,Nbase);

V1 = diag(-0.25*besselj(1,2*2*pi*Delta)*V0*ones(Nbase-1,1),1) + ...
    diag(-0.25*besselj(1,2*2*pi*Delta)*V0*ones(Nbase-1,1),-1);
V2 = diag(-0.25*besselj(2,2*2*pi*Delta)*V0*ones(Nbase-2,1),2) + ...
    diag(-0.25*besselj(2,2*2*pi*Delta)*V0*ones(Nbase-2,1),-2);

Omega = diag(repmat(omega,1,Nbase));

H = [Hbloch-Omega V1 V2; V1 Hbloch V1; V2 V1 Hbloch+Omega];
end

