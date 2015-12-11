clear all
close all
run('/Users/Eric/local/adimat-0.6.0-4971-GNU_Linux-x86_64/ADiMat_startup');

%% Run settings
global nx x dx k dt u Q nu turbmodel Cs Delta tausgs 
nx = 40;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);
k = linspace(-nx/2,nx/2 - 1,nx)';
t = 0;
et = 2;
nu = 1e-2;
dt = 1e-3;
iter = 0;

turbmodel = 2;'smagorinsky';
Cs = 0.2;
Delta = pi/20;



h = figure;
DNS = load('DNSdata/DNSkurgtad.mat');
Delta = pi/20.;
kc = 20;
%%% Initial Conditions=
%u = sin(x)';
u = real(interp1(DNS.x,DNS.ufiltsave(:,1),x)');
[nut,tausgs] = smagorinsky(u,dx);
%tausgs = interp1(DNS.x, DNS.tauSGSsave(:,1),x);
uhat = myfft(u,nx);
U2Qmap();
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
%%% Filter Info
usave = [u];
uhatsave = [uhat];
tsave = [t];
while (t < et)
  iter = iter + 1;
  advanceSolution(); 
  Q2Umap();
  t = t + dt;
  if (mod(iter,10) == 0)
    [val,DNSind] = min(abs(t - DNS.tsave));
    tsave = [tsave,t];
    usave = [usave,u];
    uhatsave = [uhatsave;uhat];
    iter2 = iter2 + 1;
    subplot(1,2,1)
    plot(x,real(u))
    hold on
    plot(DNS.x,DNS.ufiltsave(:,DNSind));
    drawnow
    hold off
    subplot(1,2,2)
    plot(x,real(tausgs))
    hold on
    %plot(DNS.x,real(DNS.tauSGSsave(:,DNSind)))
    plot(x,interp1(DNS.x,DNS.tauSGSsave(:,DNSind),x))
    hold off
    t
  end
end



