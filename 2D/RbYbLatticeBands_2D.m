clear all

NBase = 5;
V0 = 0.2;

qx = -1:0.05:1;
qy = (-1:0.05:1)';
qx = repmat(qx,length(qy),1);
qy = repmat(qy,1,length(qx));
E00 = eig(HLattice2D(0,0,V0,NBase)); E0 = E00(1);

%% Use hardcoded definition of the Hamiltionian
tic
for j = 1:length(qx)
    for l = 1:length(qy)
        E = eig(HLattice2D(qx(l,j),qy(l,j),V0,NBase));
        En0(l,j) = E(1)-E0;
        En1(l,j) = E(2)-E0;
        En2(l,j) = E(3)-E0;
    end
end
toc

%% Define k-vectors of beams and calculate Hamiltonian from the beams
%kVecs = [1 0; -1 0; 0 1];
%depths = [1 1 1];
% Find all unique 2-photon couplings
%[couplings elems] = FindCouplings(kVecs,depths);


%% Plot bands
surf(qx,qy,En0);
hold;
surf(qx,qy,En1);
hold;