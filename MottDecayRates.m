clear all

%% Unless otherwise stated, the units are as follows.
% - Momenta in recoil momenta.
% - Energies in recoil energies.
% - Distances are in lattice spacings.
% - As a result exp(i*k*x) -> exp(i*pi*k*x)

%% A bunch of constants (in SI)
hbar = 1.05457148e-34; h = hbar/(2*pi); 
mRb = 87*1.66053892e-27; 
mYb = 171*1.66053892e-27;
mu = mRb*mYb/(mRb + mYb);   
aRR = 5.31e-9;
gRR = 4*pi*(hbar^2)*aRR/mRb;
aYY = 0;
gYY = 4*pi*(hbar^2)*aYY/mYb;
aRY = -3.1e-9;
gYR = 4*pi*(hbar^2)*aRY/mu;


%% Set resolution and lattice parameters.
Lambda = 423e-9;            % Lattice wavelength (nm)
kL = 2*pi/Lambda;           % Recoil momentum
a = Lambda/2;               % Lattice constant (nm)
Er = (hbar*kL)^2/(2*mYb);   % Recoil energy
Nbase = 21;                 % Number of plane waves to expand with
nL = 150;                   % Number of lattice sites
dk = 2/nL;                  % Resolution in momentum space
dx = 0.01;                  % Resolution in physical space
V0 =  20;                   % Lattice depth

%% Diagonalize the Hamiltonian and get state vectors and band energies.
k = -1:dk:1;                 
% Discretized values of k.
En = zeros(Nbase,length(k));        % Energy(Band,k)
C = zeros(Nbase,Nbase,length(k));   % C(BrillouinZone,Band,k)

for j=1:length(k)
    H = HLattice(k(j),V0,Nbase);    % Define Hamiltonian
    [c,ex] = eig(H);                % Diagonalize
    En(:,j) = diag(ex);             % Store band energies
    C(:,:,j) = c;                   % Store state vectors
end

% Get nearest neighbor tunneling energies.
for n=1:Nbase
    t(n) = sum(En(n,:).*arrayfun(@exp,i*pi*k));
end

% This reorients the state vector to make it more usable later.
C = permute(C,[3,1,2]);             % C(k,BrllouinZone,Band)
C = C.*(2*repmat(sum(C,2)>zeros(length(k),1,Nbase),[1,Nbase,1])-1);
clear j c ex

%% Use the state vectors to get the lattice cell functions.
x = -5:dx:5;
BlochCell = zeros(length(k),length(x),Nbase);
BlochWave = BlochCell;
ex = zeros(Nbase,length(x));

% Create the plane wave matrix.
for l=1:Nbase
  for j=1:length(x)
     ex(l,j) = exp(i*pi*2*(l-(Nbase+1)/2)*x(j));    % Plane wave matrix
  end
end

% Evaluate the Bloch Cell functions using the matrix multiplication method.
for n=1:Nbase
  BlochCell(:,:,n) = C(:,:,n)*ex;                   % BlochCell(k,x,Band)
end

clear l j ex

% Get the full Bloch function by sticking on the exp(ikx) phase.
% Note that this uses element-wise multiplication, not matrix
% multiplication
ex = arrayfun(@exp,i*pi*transpose(k)*x);        % The phase matrix exp(k,x)
for n=1:Nbase
    BlochWave(:,:,n) = BlochCell(:,:,n).*ex;
end

clear n ex

%% Using the Bloch functions, calculate the Wannier functions.
Wannier = zeros(Nbase,length(x));
U = zeros(1,Nbase);

for n=1:Nbase
    Wannier(n,:) = sum(BlochWave(:,:,n),1)./(length(k));
end

% Use Wannier functions to get onsite interaction energies in the ground 
% band.
for n=1:Nbase
    U(n) = (gYY*a^3)*sum(Wannier(n,:).^4.*dx);
end
