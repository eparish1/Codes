clear all
close all


nx = 4000;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);

k = linspace(-nx/2,nx/2 - 1,nx)';
%u = sin(x)';
u = initialCondition(x,20,k);
%% Get Fourier Coefficients and wave numbers
%uhat = F*u;
uhat = myfft(u,nx);
h = figure;

Delta = 100.*dx;
Gf = heaviside(pi/Delta - abs(k));

t = 0;
et = 2;
nu = 1e-2;
%dt = nu*dx^2*20;
dt = 1e-4;
iter = 0;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
%%% Filter Info
Delta = 100.*dx;
usave = [u];
uhatsave = [uhat];
ufiltsave = [myifft(Gf.*uhat,nx)];
tauSGSsave = [myifft(computeSGS(Gf,uhat),nx)];
tsave = [t];
while (t < et)
  iter = iter + 1;
  u0 = u;
  u20 = u;
  for L=1:4
   [RHS1] = kurgTadmorRHS(u,nu,dx);
    u = u0 + dt*rk4const(L)*RHS1;
   %[RHS2] = upwindRHS(u2,nu,dx);
   % u2 = u20 + dt*rk4const(L)*RHS2;
 end
  t = t + dt;
  if (mod(iter,10) == 0)
    tsave = [tsave,t];
    usave = [usave,u];
    uhat = myfft(u,nx);
    uhatsave = [uhatsave;uhat];
    ufilt = myifft(Gf.*uhat,nx);
    ufiltsave = [ufiltsave,ufilt];
    tauSGSsave = [tauSGSsave,myifft(computeSGS(Gf,uhat),nx)];
    iter2 = iter2 + 1;
    %subplot(1,2,1)
    %plot(x,u)
    %hold on 
    %plot(x,ufilt)
    %drawnow
    %hold off
    t
  end
end

save('DNSdata/DNSkurgtad.mat','usave','uhatsave','tauSGSsave','k','x','tsave','ufiltsave')


