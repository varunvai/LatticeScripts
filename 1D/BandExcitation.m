clear all

%% Unless otherwise stated, the units are as follows.
% - Momenta in recoil momenta.
% - Energies in recoil energies.
% - Distances are in lattice spacings.
% - As a result exp(i*k*x) -> exp(i*pi*k*x)

%% A bunch of constants (in SI)
hbar = 1.05457148e-34; h = 2*pi*hbar; 
mRb = 87*1.66053892e-27; 
mYb = 171*1.66053892e-27;
mu = mRb*mYb/(mRb + mYb);   
aRR = 5.31e-9;
gRR = 4*pi*(hbar^2)*aRR/mRb;
aRY = -3.1e-9;
gYR = 4*pi*(hbar^2)*aRY/mu;

%% Set resolution and lattice parameters.
Lambda = 423.018e-9;            % Lattice wavelength (m)
kL = 2*pi/Lambda;           % Recoil momentum
a = Lambda/2;               % Lattice constant (m)
ErYb = (hbar*kL)^2/(2*mYb); % Yb Recoil energy (J)
Nbase = 41;                 % Number of plane waves to expand with
V0Yb = 25;                  % Lattice depth
theta = (0.01)*2*pi;        % Phase modulation depth

%% Diagonalize the Hamiltonian and get state vectors and band energies.
%dk=0.01;
k = -1:0.01:1;                      % Discretized values of k.
EnYb = zeros(Nbase);        % Energy(Band,k)
CYb = zeros(Nbase,Nbase);   % C(BrillouinZone,Band,k)

for j = 1:length(k)
    HYb = HLattice(k(j),V0Yb,Nbase);
    [cYb,exYb] = eig(HYb);
    H_Bloch(:,:,j) = exYb;
    
    cYbD = circshift(cYb,1);
    cYbU = circshift(cYb,-1);
    V_PM(:,:,j) = (V0Yb*theta/2)*((cYbD'*cYb-cYbU'*cYb)/2*i);
end

%% Get Rabi frequencies from ground to first band (in Hz)
for j=1:length(k)
    Omega12(j) = abs(V_PM(1,2,j));
end

%% Calculate transition matrix elements for phase modulated lattice.
%CYbd = circshift(CYb,1);
%CYbu = circshift(CYb,-1);
%V = (V0Yb*theta/2)*((CYbd'*CYb-CYbu'*CYb)/2*i);