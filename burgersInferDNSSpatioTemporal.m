clear all
close all
%run('/home/eric/local/adimat-0.6.0-4971-GNU_Linux-x86_64/ADiMat_startup');
run('/Users/Eric/local/adimat-0.6.0-4971-GNU_Linux-x86_64/ADiMat_startup');
h = figure;
%% Grid
nx = 40;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);
k = linspace(-nx/2,nx/2 - 1,nx)';

%% Initial Condition
%u = initialCondition(x,20,k);
%uhat = fftshift(fft(u))/1000.;

%% Filters
Delta = pi/20.;
kc = 20;
%Gf = exp(-k.^2*Delta^2./24.);%Gauss Filter
Gf = heaviside(pi/Delta - abs(k));%Sharp Spectral
%% Load DNS Data

DNS = load('DNSdata/DNS.mat');
GfDNS = heaviside(pi/Delta - abs(DNS.k));%Sharp Spectral
%% Initial Condition
uhat = interp1(DNS.k,GfDNS.*DNS.uhatsave(:,1),k);
u = myifft(uhat,nx);
%u = sin(x)';%interp1(DNS.x,myifft(GfDNS.*DNS.uhatsave(:,1),length(DNS.x)),x)';
%uhat = myfft(u,nx);
%u = initialCondition(x,20,k);
%uhat = myfft(u,nx);
%uhat = LES.uhatsave(:,1);
%uhat_beta = zeros(size(uhat));
t = 0;
et = 1;
nu = 1e-2;
%dt = nu*dx^2*20;
dt0 = 1e-3;
dt = 1e-3;
iter = 1;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
uhatsave(:,iter) = uhat;
tauhatsave(:,iter) = EDQNM(k,kc,uhat);
Ehatsave(:,iter) = uhat.*conj(uhat);
tsave(:,iter) = 0;
mydraw = 1;



%%% Linear Tangent Info
assimilation_window = 1%length(LES.uhatsave(1,:));
nparams = nx; %number of inferred parameters
uhat_beta = zeros(nx,nparams);
M = eye(nx*assimilation_window);
Caa = eye(nparams)*1e-3;
%Caa(1,2) = -0.1;
%Caa(2,1) = -0.1;
%Caa(2,2) = 0.2*Caa(1,1);
%Caa(2,2) = 1*Caa(1,1);
Wee = eye(nx*assimilation_window).*1./(1.e-10);
opt_its = 300;
gam = -1e-5;
beta0 = ones(size(u)).*0.5;
beta = beta0;
beta_save = beta0;

uhat_betasave = norm([uhat_beta]);
uhatsave_init = [uhat];
tsave_init = [0]
% Get initial Baseline solve

while (t < et)
  [RHS1,tauhat] = convolutionRHS_LES_beta(uhat,k,kc,nu,beta);
  Jac = admDiffFD(@convolutionRHS_LES_beta,1,uhat,k,kc,nu,beta,admOptions('i',1,'d',1));
  A = eye(nx)./dt - Jac;
  uhat = A\RHS1 + uhat;
  uhatsave_init = [uhatsave_init,uhat];
  t = t + dt
  tsave_init = [tsave_init,t];
end
t = 0;
uhat = interp1(DNS.k,GfDNS.*DNS.uhatsave(:,1),k);
while (t < et)
  iter = iter + 1;
  uhat0 = uhat;
  uhat_beta = zeros(size(uhat_beta));
  uhat_beta0 = uhat_beta;
  for OPT_ITER=1:opt_its
    [RHS1,tauhat] = convolutionRHS_LES_beta(uhat0,k,kc,nu,beta);
    Jac = admDiffFD(@convolutionRHS_LES_beta,1,uhat0,k,kc,nu,beta,admOptions('i',1,'d',1));
    dGdbeta =  admDiffFD(@convolutionRHS_LES_beta,1,uhat0,k,kc,nu,beta,admOptions('i',5,'d',1));
    uhat_beta = uhat_beta0 + dt*((Jac*uhat_beta) - dGdbeta);
    A = eye(nx)./dt - Jac;
    uhat = A\RHS1 + uhat0;
    [val,DNSind] = min(abs(t - DNS.tsave));
    obs_error =  interp1(DNS.k,DNS.uhatsave(:,DNSind),k) - uhat;
    beta = beta - real(gam*( -Caa * (M * uhat_beta)'*Wee*obs_error));
  end
  uhat_betasave = [uhat_betasave;norm(uhat_beta(:,1))];
  beta_save = [beta_save,beta];
  uhatsave = [uhatsave;uhat];
  t = t + dt;
  if (mod(iter,10) == 0)
    [val,DNSind] = min(abs(t - DNS.tsave));
    [val,iniInd] = min(abs(t - tsave_init));
    iter2 = iter2 + 1;
    tauhatsave(:,iter2) = tauhat(:,1);
    Ehatsave(:,iter2) = uhat.*conj(uhat);
    tsave(iter2) = t;
    if (mydraw == 1)
      subplot(1,2,1)
      plot(abs(k),abs(uhat),'color','blue')
      hold on
      plot(abs(k),abs(uhatsave_init(:,iniInd)),'color','red')
      %loglog(abs(k),abs(tauhat))
      plot(abs(DNS.k),abs(GfDNS.*DNS.uhatsave(:,DNSind)))
      xlim([0,20])
      hold off
      %loglog(abs(LES.k),abs(LES.tauhatsave(:,LESind)),'o','color','red')
      %%plot(ones(2,1)*kc,linspace(1e-10,1,2),'color','black','LineWidth',3)
      %hold off
      %xlabel('$k$','Interpreter','LaTex','FontSize',20)
      %ylabel('$|\tau_{sgs}|$','Interpreter','LaTex','FontSize',20)
      %legend('DNS','Filtered DNS','Filter Cutoff')
      
      subplot(1,2,2)
      plot(x,real(myifft(uhat,nx)),'linewidth',2)
      %loglog(abs(k),uhat.*conj(uhat))
      hold on
      plot(DNS.x,real(myifft(GfDNS.*DNS.uhatsave(:,DNSind),length(DNS.k))))
      %loglog(abs(LES.k),LES.Ehatsave(:,LESind),'color','red')
      %plot(ones(2,1)*kc,linspace(1e-10,1,2),'o','color','black','LineWidth',3)
      xlabel('$k$','Interpreter','LaTex','FontSize',20)
      ylabel('$\hat{uu^{\ast}}$','Interpreter','LaTex','FontSize',20)
      legend('LES','Filtered DNS')
      %ylim([1e-10,1e8])
      %xlim([1,5.*max(k)])
      hold off
      drawnow
    end
    t
    beta
  end
end

DNSobs = [];
for i = 1:1001
  DNSobs = [DNSobs;DNS.uhatsave(end/2+1-kc:end/2+kc,i)];
end
u2 = reshape(uhatsave_init,[40040,1]);
u3 = reshape(uhatsave,[40040,1]);

%close all
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0,'LESFigures/edq_tau000.png',1)
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0.1,'LESFigures/edq_tau010.png',1)
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0.2,'LESFigures/edq_tau020.png',1)
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0.3,'LESFigures/edq_tau030.png',1)
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0.4,'LESFigures/edq_tau040.png',1)
%plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,0.5,'LESFigures/edq_tau050.png',1)

%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0,'LESFigures/edq_E000.png',1)
%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0.1,'LESFigures/edq_E010.png',1)
%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0.2,'LESFigures/edq_E020.png',1)
%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0.3,'LESFigures/edq_E030.png',1)
%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0.4,'LESFigures/edq_E040.png',1)
%plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,0.5,'LESFigures/edq_E050.png',1)
