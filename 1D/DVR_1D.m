clear all

h = 6.62607004e-34;
hbar = h/(2*pi);
m = 171*1.66054e-27;

Ngrid = 2048;
N = Ngrid+2;
Lgrid = 4.23018e-6;
x = linspace(-Lgrid/2,Lgrid/2,Ngrid); 
dx = x(2)-x(1);
L = Lgrid+2*dx;

lambda = 423.018e-9;
kL = 2*pi/(423.018e-9);
Er = (hbar*kL)^2/(2*m);
V0 = 30*Er;
a = 0.05;
pot = (V0/2)*(1-cos(2*kL.*x)) + a*Er*x/lambda;
%pot = 0.5*m*(2*pi*100)^2*x.^2;
V = diag(pot)/h;
theta = (0.01)*2*pi;

T = zeros(Ngrid,Ngrid);
j = repmat((1:Ngrid)',1,Ngrid);
l = repmat((1:Ngrid),Ngrid,1);
T = ((hbar*pi)^2/(4*m*L^2))*((-1).^(j-l)).*(sin(pi*(j-l)/(2*N)).^-2 ...
    + sin(pi*(j+l)/(2*N)).^-2);
T(logical(eye(size(T))))=0;
T = T + diag((hbar*pi)^2/(4*m*L^2)*((2*N^2+1)/3 - 1./(sin(pi*(1:Ngrid)/N).^2)));
T = T/h;

H = T+V;
[V,En] = eig(H);

    