clear all

Nbase = 5;
V0 = 7/9;
Omega = 2.2;
Delta = 0.01;

qx = -1:0.05:1;
qy = (-1:0.05:1)';
qx = repmat(qx,length(qy),1);
qy = repmat(qy,1,length(qx));
E00 = eig(HFloquet2D(0,0,V0,Delta,Omega,Nbase)); E0 = E00(1);
En0 = qy; En1 = qy;

tic
for j = 1:length(qx)
    for l = 1:length(qy)
        H = HFloquet2D(qx(l,j),qy(l,j),V0,Delta,Omega,Nbase);
        E = eig(H); 
        En0(l,j) = E(1)-E0; 
        En1(l,j) = E(2)-E0;
    end
end
toc
figure(2)
surf(qx,qy,En0);
hold
surf(qx,qy,En1);
hold
