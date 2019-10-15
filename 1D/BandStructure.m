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
Lambda = 423e-9;            % Lattice wavelength (m)
kL = 2*pi/Lambda;           % Recoil momentum
a = Lambda/2;               % Lattice constant (m)
ErRb = (hbar*kL)^2/(2*mRb); % Rb Recoil energy (J)
ErYb = (hbar*kL)^2/(2*mYb); % Yb Recoil energy (J)
Nbase = 15;                 % Number of plane waves to expand with
nL = 100;                   % Number of lattice sites
dk = 2/nL;                  % Resolution in momentum space
dx = 0.01;                  % Resolution in physical space
V0Rb = 10;                     % Lattice depth
V0Yb = 20;


%% Diagonalize the Hamiltonian and get state vectors and band energies.
dk=0.01;
k = -1:dk:1;                        % Discretized values of k.
EnRb = zeros(Nbase,length(k));        % Energy(Band,k)
CRb = zeros(Nbase,Nbase,length(k));   % C(BrillouinZone,Band,k)
EnYb = EnRb; CYb = CRb;

for j=1:length(k)
    HRb = HLattice(k(j),V0Rb,Nbase);      % Rb Hamiltonian
    HYb = HLattice(k(j),V0Yb,Nbase);      % Yb Hamiltonian
    [cRb,exRb] = eig(HRb);                % Diagonalize
    [cYb,exYb] = eig(HYb);
    EnRb(:,j) = diag(exRb);               % Store band energies
    EnYb(:,j) = diag(exYb);
    CRb(:,:,j) = cRb;                     % Store state vectors
    CYb(:,:,j) = cYb; 
end
EnRb = EnRb*ErRb/h; EnYb = EnYb*ErYb/h;

%% BEC spectrum
%rho = 10^20; kdense = -3:0.01:3;
%EbecArr = sqrt(((gRR*rho/mRb)*(hbar*kdense*kL).^2) + ((hbar*kdense*kL).^2/(2*mRb)).^2)/Er;


%E0 = polyfit(k,En(1,:)-En(1,ceil(length(k)/2)),18);
%E1 = polyfit(k,En(2,:)-En(1,ceil(length(k)/2)),30);
%Ebec = polyfit(kdense,EbecArr,50);

figure(1);
plot(k,EnRb(1,:)-EnRb(1,ceil(length(k)/2)),k,EnRb(2,:)-EnRb(1,ceil(length(k)/2))...
,k,EnYb(1,:)-EnYb(1,ceil(length(k)/2)),k,EnYb(2,:)-EnYb(1,ceil(length(k)/2)));

figure(2);
plot(k,EnRb(1,:)-EnRb(1,ceil(length(k)/2)),k,EnYb(1,:)-EnYb(1,ceil(length(k)/2)));