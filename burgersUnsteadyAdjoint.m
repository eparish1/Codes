clear all
close all
run('/Users/Eric/local/adimat-0.6.0-4971-GNU_Linux-x86_64/ADiMat_startup');
DNS = load('Case2/DNS/DNS.mat');
global nx x dx k dt u uhat Q nu turbmodel Cs Delta tausgs tauhatsgs convective_flux solution_domain M t Gf dQdBeta beta 
%% ============ SETTINGS ==========================
solution_domain = 2; %1 for spatial and 2 for frequency
convective_flux = 2; %1 for upwind and 2 for MUSCL
turbmodel = 5;     %0 for DNS, 1 for smag, 2 for lag, 3 for exact
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
et = 2e-3;
dt = 1;
M = linspace(-nx,nx-1,2*nx);
for i = 1:2*nx
  if (M(i) <= nx/2 && M(i) >= -nx/2)
    M(i) = 0; 
  end
end
M = M';
%%==============================================


%% ============= Initial Conditions  =================
GfDNS = heaviside(pi/Delta - abs(DNS.k));
Gf = heaviside(pi/Delta - k);
%u = sin(x)';
[dum,j_min] = min(abs(DNS.k - -kc));
[dum,j_max] = min(abs(DNS.k - (kc-1)));
[dum,i_min] = min(abs(k - -kc));
[dum,i_max] = min(abs(k - (kc-1)));
uhat = zeros(nx,1);
%u = initialCondition(x,kc,k);
%uhat = myfft(u,nx);
uhat(i_min:i_max,1) =  DNS.uhatsave(j_min:j_max,1);
u = myifft(uhat,nx);
%u = real(interp1(DNS.x,myifft(GfDNS.*DNS.uhatsave(:,1),length(DNS.k)),x)');
[nut,tausgs] = smagorinsky(u,dx);
%tausgs = interp1(DNS.x, DNS.tauSGSsave(:,1),x);
%uhat = myfft(u,nx);
U2Qmap();
QIC = Q;
nu = 1e-2;
uIC = u;
uhatIC = uhat;
%%===================================================


%%=============== Inversion setup ==================
assimilation_window = 1;%length(DNS.uhatsave(1,:));
nparams = nx; %number of inferred parameters
dQdBeta = zeros(length(Q),nparams); %gradient of Q w.r.p to beta
H = eye(nx*assimilation_window); %maping of Q to obs
Cbb = eye(nparams); %covariance of the prior
Wee = eye(nx*assimilation_window).*1./(1.e-10); %inverse covariance of obs
opt_its =50; % Number of optimization runs
gam = -1e-8; % step size for steepest decent
beta0 = ones(nparams,1); %initial beta
beta = beta0; %
beta_save = beta0;
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
error_save = [];
while (t < et)
  iter = iter + 1;
  Q0 = Q;
  uhat0 = uhat;
  dQdBeta = zeros(length(Q),nparams); %gradient of Q w.r.p to beta
  beta_tmp = beta;
  for OPT_ITER = 1:opt_its
    Q = Q0;
    uhat = uhat0;
    Q = advanceSolutionUnsteadyAdjoint(Q,beta);
    %% Solve adjoint
    [val,DNSind] = min(abs(t +dt  - DNS.tsave));
    d = interp1(DNS.k,DNS.uhatsave(:,DNSind),k); 
    dhoRdhoQ = admDiffFD(@advanceSolutionUnsteadyAdjoint,1,Q,beta,admOptions('i',1,'d',1));
    dhoRdhoBeta = admDiffFD(@advanceSolutionUnsteadyAdjoint,1,Q,beta,admOptions('i',2,'d',1));
    dhoFdhoU = admDiffFD(@objectiveFunction,1,uhat,d,beta,beta,admOptions('i',1,'d',1));
    dhoFdhoBeta = admDiffFD(@objectiveFunction,1,uhat,d,beta,beta,admOptions('i',2,'d',1));
    psiT = dhoRdhoQ.'\(-dhoFdhoU');
    dJdBeta = dhoFdhoBeta + psiT.'*dhoRdhoBeta;
    beta = beta + 1e-5*dJdBeta';
    beta_tmp = [beta_tmp,beta];
    Q2Umap();
%    [val,DNSind] = min(abs(t - DNS.tsave));
%    obs_error = interp1(DNS.k,DNS.uhatsave(:,DNSind),k)  - uhat;
%    beta = beta - real(gam.*( - Cbb * (H * dQdBeta)'*Wee*obs_error));
%    beta(nx/2+2:end) = real(flipud(beta(2:nx/2))) - imag(flipud(beta(2:nx/2)));
  end
  t = t + dt;
  if (mod(iter,save_freq) == 0)
    [val,DNSind] = min(abs(t - DNS.tsave));
    tsave = [tsave,t];
    usave = [usave,u];
    uhatsave = [uhatsave,uhat];
    tauhatsave = [tauhatsave,tauhatsgs];
    beta_save = [beta_save beta];
    iter2 = iter2 + 1;
    if (live_plot == 1)
      if (solution_domain == 1)
        uhat = myfft(u,nx);
      end
      if (solution_domain == 2)
       u = myifft(uhat,nx);
      end
      subplot(1,2,1)
      loglog(k,abs(uhat),'o')
      hold on
      loglog(DNS.k,abs(GfDNS.*DNS.uhatsave(:,DNSind)));
      xlim([1,50])
      drawnow
      hold off
      subplot(1,2,2)
      plot(x,myifft(tauhatsgs,nx))
      hold on
      plot(x,interp1(DNS.x,myifft(DNS.tauhatsave(:,DNSind),length(DNS.x)),x))
      %xlim([min(k),max(k)])
      hold off
      hold off
    end
    t
  end
end
%  obs_error = d  - reshape(uhatsave,[nx*assimilation_window,1]);
%  beta = beta - real(gam.*( beta - beta0 - Cbb * (H * dQdBetasave)'*Wee*obs_error));
%  %beta(nx/2+2:end) = flipud(beta(2:nx/2));
%  beta_save = [beta_save beta]
%  error_save = [error_save,norm(obs_error)]


%save('Case2/Inference/tmodel_betafullxt/Sol.mat','usave','uhatsave','tauhatsave','k','x','tsave','kc','beta_save','error_save');

