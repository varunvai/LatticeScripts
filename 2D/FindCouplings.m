function [couplings elems] = FindCouplings(kVecs,depths)
%Given a list of k-vectors (kVecs) and single beam trap depths (depths),
%find all possible 2-photon couplings kVecs(j)-kVecs(i) between pairs,
%along with associated elements (elems) of the Hamiltonian

%% Calculate all possible couplings in [M;N] basis
% For NFields beams, there are NCouplings = N!/[(N-2)!2!] couplings
kVecs = [1 0; -1 0; 0 1];
NkVecs = size(kVecs,1);
NCouplings = nchoosek(size(kVecs,1),2);
couplings = zeros(NCouplings,2);

elems = zeros(NCouplings,1);

% Find all unique couplings k2-k1 between kVecs.
n=1; ind1 = zeros(NkVecs,1); ind2 = zeros(NkVecs,1);
for j = 1:(NkVecs-1)
    for l = (j+1):NkVecs
        couplings(n,:) = kVecs(l,:)-kVecs(j,:);
        elems(n) = sqrt(depths(l))*sqrt(depths(j));
        ind2(n) = l; ind1(n) = j;
        l = l+1; n = n+1;
    end
end

end

