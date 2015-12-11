clear all
close all
%DNS = load('DNSdata/DNS.mat');
DNS = load('Case2/DNS/DNS.mat');
%% ============ SETTINGS ==========================
Cs = 0.2;            % Smagorinsky constant
kc = 200;
Delta = pi/kc;

%% ========== MESH and Time Setup =====================
nx = 400;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);
xf = linspace(0,2.*pi - dx,nx*2);
xf2 = linspace(0,2.*pi - dx,nx*2+1);
k = linspace(-nx/2,nx/2 - 1,nx)';
t = 0;
et = 2;
dt = 1e-4;
M = linspace(-nx,nx-1,2*nx);
for i = 1:2*nx
  if (M(i) <= nx/2 && M(i) >= -nx/2)
    M(i) = 0; 
  end
end
M = M';

kf = linspace(-nx,nx-1,2*nx);
Gf = heaviside(pi/Delta - abs(k + 0.1) );
Gff = heaviside(pi/Delta - abs(kf + 0.1))';
Gff(3*nx/2+1) = 1.;
%%==============================================


%% ============= initial conditions  =================
gfDNS = heaviside(pi/Delta - abs(DNS.k));
%u = sin(x)';
[dum,j_min] = min(abs(DNS.k - -kc));
[dum,j_max] = min(abs(DNS.k - (kc-1)));
[dum,i_min] = min(abs(k - -kc));
[dum,i_max] = min(abs(k - (kc-1)));
uhat = zeros(nx,1);
u = initialCondition(x,kc,k);
uhat = myfft(u,nx);
%uhat(i_min:i_max,1) =  DNS.uhatsave(j_min:j_max,1);
%u = myifft(uhat,nx);

t = 1.

xf = linspace(0,L - dx,2*nx)';
uf = interp1(x,u,xf);
uhat2 = zeros(2*nx,1);
uhat2(nx/2+1:nx/2+nx,1) = uhat;
uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
tmp = -t.*1j.*M./2.*myconv(uhat2,uhat2,length(uhat2));
tauhat1 = myconv(uhat,tmp,length(uhat));
tau1 = myifft(tauhat1,nx);

L = 0.5*Gff.*myconv(uhat2,uhat2,2*nx) - 0.5*myconv(Gff.*uhat2,Gff.*uhat2,2*nx);
Lr = myifft(L,nx);
tautest = uf.*diff1(Lr,dx/2);

L2 = interp1(xf,Lr,x);
Lb = myfft(L2,nx);
tauhat2 = myconv(uhat2,-1j.*kf'.*L);
tauhat2 = tauhat2(nx/2+1 : 3*nx/2);
%tau2 = u.*diff1( interp1(xf,myifft(L,nx),x)' ,dx);
