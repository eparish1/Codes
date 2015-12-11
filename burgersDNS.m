clear all
close all


nx = 4000;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);

k = linspace(-nx/2,nx/2 - 1,nx)';
u = initialCondition(x,20,k);
%% Get Fourier Coefficients and wave numbers
%uhat = F*u;
uhat = myfft(u,nx);
h = figure;

t = 0;
et = 2.;
nu = 1e-2;
%dt = nu*dx^2*20;
dt0 = 1e-5;
dt = 1e-5;
iter = 0;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
%%% Filter Info
Delta = 100.*dx;
%Gauss Filter
%Gf = exp(-k.^2*Delta^2./24.);
%Sharp Spectral
Gf = heaviside(pi/Delta - abs(k));
uhatsave(:,1) = uhat;
tauSGSsave(:,1) =  0.5*(Gf.*myconv(uhat,uhat,nx) - myconv(Gf.*uhat,Gf.*uhat,nx));
uhatsave(:,iter2) = uhat;
Ehatsave(:,iter2) = 0.5*uhat.*conj(uhat);
ufiltsave(:,1) = Gf.*u;
TEsave = sum(uhat.*conj(uhat))*nx;

while (t < et)  
  iter = iter + 1;
  uhat0 = uhat;
  Khat0 = 0.5*myconv(uhat0,uhat0,nx);
  tauSGS0 = 0.5*(Gf.*myconv(uhat,uhat,nx) - myconv(Gf.*uhat,Gf.*uhat,nx));
  for L=1:4
    [RHS1,c] = convolutionRHS(uhat,k,nu);
    uhat = uhat0 + dt*rk4const(L)*(RHS1);
  end
  t = t + dt;
    if (mod(iter+1,100) == 0)
      %dt = 1e-7;
    end
    if (mod(iter,100) == 0)
    iter2 = iter2 + 1;
    uhatsave(:,iter2) = uhat;
    Ehatsave(:,iter2) = 0.5*uhat.*conj(uhat);
    tsave(:,iter2) = t;
    tsave(iter2) = t;
    ufiltsave(:,iter2) = [ufiltsave,myifft(Gf.*uhat,nx)];
    %subplot(1,2,1)
    %plot(x,real(ifft(fftshift(uhat*1000))))
    %plot(x,real(G*uhat))
    %hold on
    %plot(x,real(G*uhat_f))
    %xlabel('$x$','Interpreter','LaTex','FontSize',20)
    %ylabel('$u$','Interpreter','LaTex','FontSize',20)
    %legend('DNS','Filtered DNS')
    %hold off
    %subplot(1,2,2)
    %loglog(abs(k),abs(uhat))  
    %semilogy(k,abs(c))
    %loglog(abs(k),uhat.*conj(uhat))
    %loglog(abs(k),abs(k).^(-5./3.)*3.)
    %semilogy(k,abs(c).*exp(-k.^2*(50.*dx).^2./24.))
    %xlabel('$k$','Interpreter','LaTex','FontSize',20)
    %ylabel('$\hat{u^2}$','Interpreter','LaTex','FontSize',20)
    %legend('DNS','Filtered DNS')
    %ylim([1e-40,4])
    %xlim([1,max(k)])
    %set(h,'Position',[500,500,800,300]);
    %hold off
    %drawnow
    t

    %% LES Budget
    ufilt = Gf.*uhat0;
    tauSGS = 0.5*(Gf.*myconv(uhat,uhat,nx) - myconv(Gf.*uhat,Gf.*uhat,nx));
    tauSGSsave(:,iter2) = tauSGS;
    ddt_ufilt = 1./dt*(Gf.*uhat - Gf.*uhat0);
    u_filt_budget(iter2) =norm( ddt_ufilt + 1j.*k./2.*myconv(Gf.*uhat0,Gf.*uhat0,nx) + nu*k.^2.*Gf.*uhat0 + 1j.*k.*tauSGS0 );
    %[RHS1,c] = convolutionRHS(uhat,k,nu);
    u_budget(iter2) = norm(1./dt*(uhat - uhat0) - RHS1 );
    %%% SGS Equation Budget
    phi = 2./6.*(1j.*k.*myconv(ufilt,Gf.*myconv(uhat0,uhat0,nx),nx) - ...
         1j.*k.*Gf.*myconv(uhat0,myconv(uhat0,uhat0,nx),nx)) + ...
         nu.*myconv(1j.*k.*ufilt,1j.*k.*ufilt,nx) - ...
         nu.*Gf.*(myconv(1j.*k.*uhat0,1j.*k.*uhat0,nx));
   
    T1Save(:,iter2) = 2./3.*1j.*k.*myconv(ufilt,tauSGS0,nx) ;
    T2Save(:,iter2) = k.^2*nu.*tauSGS0;
    T3Save(:,iter2) = myconv(ufilt,1j.*k.*tauSGS0,nx);
    T4Save(:,iter2) = phi;
    tauSGS_budget(iter2) = norm( 1./dt*(tauSGS - tauSGS0) + T1Save(:,iter2) + ...
                           T2Save(:,iter2) - T3Save(:,iter2) - phi);
    dt = dt0;
    Khat = 0.5*myconv(uhat,uhat,nx);
    TKEbudget(iter2) = norm( 1./dt*(Khat - Khat0) + 2./3.*1j.*k.*myconv(uhat0,Khat0,nx) + nu.*k.^2.*Khat0 + ...
                       nu*myconv(1j.*k.*uhat0,1j.*k.*uhat0,nx));
    knorm(iter2) = norm(Khat);
    unorm(iter2) = norm(uhat);
    end
end
save('DNSdata/DNS.mat','uhatsave','tauSGSsave','k','x','tsave','ufiltsave')

%{
for i=1:iter2
    T1norm(i) = norm(abs(T1Save(:,i)));
    T2norm(i) = norm(abs(T2Save(:,i)));
    T3norm(i) = norm(abs(T3Save(:,i)));
    T4norm(i) = norm(abs(T4Save(:,i)));
end
%%%%%%%%%%%%%%%
[t0,t0ind] = min(abs(tsave - 0.00)) 
t0ind = t0ind + 1;
[t1,t1ind] = min(abs(tsave - 0.05))
[t2,t2ind] = min(abs(tsave - 0.1))
[t3,t3ind] = min(abs(tsave - 0.15))
close all
Emax = max(Ehatsave(nx/2+2:end,:))
Emin = [1e-30;1e-30];%min(Ehatsave(nx/2+2:end,:))*10


figure(1)
plot(x,real(myifft(uhatsave(:,t0ind),nx)),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$x$','Interpreter','LaTex','FontSize',30)
ylabel('$u$','Interpreter','LaTex','FontSize',30)
saveas(gcf,'DNSFigures/t0u.png')


figure(2)
loglog(k,Ehatsave(:,t0ind),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$k$','Interpreter','LaTex','FontSize',30)
ylabel('$\hat{u^2}$','Interpreter','LaTex','FontSize',30)
ylim([Emin(2)/2,Emax(2)*2])
xlim([1,max(k)])
saveas(gcf,'DNSFigures/t0E.png')

close all 
figure(1)
plot(x,real(myifft(uhatsave(:,t1ind),nx)),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$x$','Interpreter','LaTex','FontSize',30)
ylabel('$u$','Interpreter','LaTex','FontSize',30)
saveas(gcf,'DNSFigures/t1u.png')
figure(2)
loglog(k,Ehatsave(:,t1ind),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$k$','Interpreter','LaTex','FontSize',30)
ylabel('$\hat{u^2}$','Interpreter','LaTex','FontSize',30)
ylim([Emin(2)/2,Emax(2)*2])
xlim([1,max(k)])
saveas(gcf,'DNSFigures/t1E.png')

figure(3)
plot(x,real(myifft(uhatsave(:,t2ind),nx)),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$x$','Interpreter','LaTex','FontSize',30)
ylabel('$u$','Interpreter','LaTex','FontSize',30)
saveas(gcf,'DNSFigures/t2u.png')

figure(4)
loglog(k,Ehatsave(:,t2ind),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$k$','Interpreter','LaTex','FontSize',30)
ylabel('$\hat{u^2}$','Interpreter','LaTex','FontSize',30)
ylim([Emin(2)/2,Emax(2)*2])
xlim([1,max(k)])
saveas(gcf,'DNSFigures/t2E.png')

figure(5)
plot(x,real(myifft(uhatsave(:,t3ind),nx)),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$x$','Interpreter','LaTex','FontSize',30)
ylabel('$u$','Interpreter','LaTex','FontSize',30)
saveas(gcf,'DNSFigures/t3u.png')

figure(6)
loglog(k,Ehatsave(:,t3ind),'color','blue','linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$k$','Interpreter','LaTex','FontSize',30)
ylabel('$\hat{u^2}$','Interpreter','LaTex','FontSize',30)
ylim([Emin(2)/2,Emax(2)*2])
xlim([1,max(k)])
saveas(gcf,'DNSFigures/t3E.png')


figure(7)
semilogy(tsave,u_budget,'o','color','blue')
hold on
semilogy(tsave,u_filt_budget,'linewidth',2,'color','red')
semilogy(tsave,tauSGS_budget,'linewidth',2,'color','black')
hold off
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$t$','Interpreter','LaTex','FontSize',30)
ylabel('$Budget$','Interpreter','LaTex','FontSize',30)
legend('DNS','LES','SGS')
saveas(gcf,'DNSFigures/budgets.png')

figure(8)
semilogy(tsave,T1norm,'linewidth',2,'color','blue')
hold on
semilogy(tsave,T2norm,'linewidth',2,'color','red')
semilogy(tsave,T3norm,'linewidth',2,'color','green')
semilogy(tsave,T4norm,'linewidth',2,'color','black')
hold off
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$t$','Interpreter','LaTex','FontSize',30)
ylabel('$SGS$ $Budget$','Interpreter','LaTex','FontSize',30)
legend('T_1','T_2','T_3','T_4')
saveas(gcf,'DNSFigures/SGSbudgets.png')

figure(9)
plot(tsave,T4norm./(T1norm + T2norm + T3norm + T4norm)*100,'linewidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',20)
xlabel('$t$','Interpreter','LaTex','FontSize',30)
ylabel('$\% of Budget$','Interpreter','LaTex','FontSize',30)
saveas(gcf,'DNSFigures/SGSbudgetsPercent.png')

%}

