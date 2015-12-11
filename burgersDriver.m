clear all
close all
run('/Users/Eric/local/adimat-0.6.0-4971-GNU_Linux-x86_64/ADiMat_startup');
%DNS = load('DNSdata/DNS.mat');
DNS = load('Case2/DNS/DNS.mat');
global nx x dx k dt u uhat Q nu turbmodel Cs Delta tausgs tauhatsgs convective_flux solution_domain ...
       M t Gf w0hat w1hat w2hat
%% ============ SETTINGS ==========================
solution_domain = 2; %1 for spatial and 2 for frequency
convective_flux = 2; %1 for upwind and 2 for MUSCL
turbmodel = 14;     %0 for DNS, 1 for smag, 2 for lag, 3 for exact
Cs = 0.2;            % Smagorinsky constant
save_freq = 1;     % Frequency to save solution to usave
live_plot = 1;       % Live plotting of solution
%%% Filter Info
Delta = pi/20;
kc = 20;
Delta = pi/20.;

%% ========== MESH and Time Setup =====================
nx = 40;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);
k = linspace(-nx/2,nx/2 - 1,nx)';
t = 0;
et = 2;
dt = 1e-3;
M = linspace(-nx,nx-1,2*nx);
for i = 1:2*nx
  if (M(i) <= nx/2 && M(i) >= -nx/2)
    M(i) = 0; 
  end
end
M = M';
Gf = heaviside(pi/Delta - abs(k + 0.1) );
%%==============================================


%% ============= Initial Conditions  =================
GfDNS = heaviside(pi/Delta - abs(DNS.k));
%u = sin(x)';
[dum,j_min] = min(abs(DNS.k - -kc));
[dum,j_max] = min(abs(DNS.k - (kc-1)));
[dum,i_min] = min(abs(k - -kc));
[dum,i_max] = min(abs(k - (kc-1)));
uhat = zeros(nx,1);
w0hat = zeros(nx,1);
w1hat = zeros(nx,1);
w2hat = zeros(nx,1);

%u = initialCondition(x,kc,k);
%u = sin(x)';
%uhat = myfft(u,nx);
uhat(i_min:i_max,1) =  DNS.uhatsave(j_min:j_max,1);
u = myifft(uhat,nx);
%u = real(interp1(DNS.x,myifft(GfDNS.*DNS.uhatsave(:,1),length(DNS.k)),x)');
[nut,tausgs] = smagorinsky(u,dx);
tauhatsgs = zeros(nx,1);
%tausgs = interp1(DNS.x, DNS.tauSGSsave(:,1),x);
%uhat = myfft(u,nx);
U2Qmap();
nu = 1e-2;
%%===================================================

%% other initializations
h = figure;
iter = 0;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
usave = [u];
uhatsave = [uhat];
tsave = [t];
tauhatsave = myfft(tausgs,nx);


%%%% =============== 

while (t < et)
  iter = iter + 1;
  advanceSolution(); 
  Q2Umap();
  t = t + dt;
  if (mod(iter,save_freq) == 0)
    [val,DNSind] = min(abs(t - DNS.tsave));
    tsave = [tsave,t];
    usave = [usave,u];
    uhatsave = [uhatsave,uhat];
    tauhatsave = [tauhatsave,tauhatsgs];
    iter2 = iter2 + 1;
    if (live_plot == 1)
      if (solution_domain == 1)
        uhat = myfft(u,nx);
      end
      if (solution_domain == 2)
       u = myifft(uhat,nx);
      end
      subplot(1,2,1)
      plot(x,real(u))
      hold on
      plot(DNS.x,real(myifft(DNS.uhatsave(:,DNSind),length(DNS.k))));
      drawnow
      hold off
      subplot(1,2,2) 
      plot(k,abs(tauhatsgs))
%      xlim([-100,100])
%      plot(x,real(myifft(tauhatsgs,nx)))
      hold on
%      plot(x,real(interp1(DNS.x,myifft(DNS.tauhatsave(:,DNSind),length(DNS.x)),x)))
%      legend('Model','DNS')
      val = 0.5*(GfDNS.*myconv(DNS.uhatsave(:,DNSind),DNS.uhatsave(:,DNSind),nx)  - myconv(GfDNS.*DNS.uhatsave(:,DNSind),GfDNS.*DNS.uhatsave(:,DNSind),nx));
      plot(k,abs(interp1(DNS.k,val,k)))
%      plot(x,real(interp1(DNS.x,myifft(DNS.tauhatsave(:,DNSind),length(DNS.k)),x)))
      %xlim([min(k),max(k)])
      hold off
      hold off
    end
    t
  end
end

%save('Case3/FM1/LES.mat','usave','uhatsave','tauhatsave','k','x','tsave','kc');



