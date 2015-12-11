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
tauhat_DNS = computeSGS(GfDNS,DNS.uhatsave(:,1));
tauhat = interp1(DNS.k,tauhat_DNS,k);

%% Initial Condition
%uhat = interp1(DNS.k,DNS.uhatsave(:,1),k);
%u = myifft(uhat,nx);
%u = sin(x)';
u = interp1(DNS.x,myifft(GfDNS.*DNS.uhatsave(:,1),length(DNS.x)),x)';
uhat = myfft(u,nx);
%u = initialCondition(x,20,k);
%uhat = myfft(u,nx);
t = 0;
et = 0.1;
nu = 1e-2;
%dt = nu*dx^2*20;
dt0 = 1e-4;
dt = 1e-4;
iter = 1;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
uhatsave(:,iter) = uhat;
tauhatsave(:,iter) = EDQNM(k,kc,uhat);
Ehatsave(:,iter) = uhat.*conj(uhat);
tsave(:,iter) = t;
mydraw = 1;
while (t < et)
  
  iter = iter + 1;
  uhat0 = uhat;
  for L=1:4
    [RHS1,tauhat] = convolutionRHS_LES(uhat,tauhat,k,kc,nu);
  %  Jac = admDiffFD(@convolutionRHS_LES,1,uhat,tauhat,k,kc,nu,admOptions('i',1,'d',1));
    uhat = uhat0 + dt*rk4const(L)*(RHS1);
  end
  %[RHS1,tauhat] = convolutionRHS_LES(uhat,tauhat,k,kc,nu);
  %Jac = admDiffFD(@convolutionRHS_LES,1,uhat,tauhat,k,kc,nu,admOptions('i',1,'d',1));
  %A = eye(nx)./dt - Jac;
  %uhat = A\RHS1 + uhat;
  t = t + dt;
  tauhatsave(:,iter) = tauhat;
  uhatsave(:,iter) = uhat;
  Ehatsave(:,iter) = uhat.*conj(uhat);
  tsave(iter) = t;
  if (mod(iter,100) == 0)
    %% DNS Index:
    [val,DNSind] = min(abs(t - DNS.tsave));
    
    iter2 = iter2 + 1;
    if (mydraw == 1)
      subplot(1,2,1)
      plot(x,myifft(uhat,nx),'color','blue')
      hold on
      plot(DNS.x,myifft(GfDNS.*DNS.uhatsave(:,DNSind),length(DNS.x)),'color','red')
      hold off
      xlabel('$k$','Interpreter','LaTex','FontSize',20)
      ylabel('$|\tau_{sgs}|$','Interpreter','LaTex','FontSize',20)
      legend('LES','Filtered DNS','Filter Cutoff')
      
      subplot(1,2,2)
      %plot(x,myifft(uhat,nx))
      loglog(abs(k),uhat.*conj(uhat))
      hold on
      loglog(abs(DNS.k),GfDNS.*DNS.uhatsave(:,DNSind).*conj(GfDNS.*DNS.uhatsave(:,DNSind)),'color','red')
      plot(ones(2,1)*kc,linspace(1e-10,1,2),'color','black','LineWidth',3)
      xlabel('$k$','Interpreter','LaTex','FontSize',20)
      ylabel('$\hat{uu^{\ast}}$','Interpreter','LaTex','FontSize',20)
      legend('LES','Filtered DNS','Filter Cutoff')
      %ylim([1e-10,1e8])
      xlim([1,5.*max(k)])
      set(h,'Position',[500,500,800,300]);
      hold off
      drawnow
    end
    t
  end
end
save('LESdata/LES.mat','uhatsave','tauhatsave','k','x','tsave','Ehatsave')
close all
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
