function [BlochPsi] = Bloch(x,k,V0,n,Nbase)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

H = HLattice(k,V0,Nbase);
[C,E]=eig(H);
c=C(:,1);

BlochPsi = 
end

