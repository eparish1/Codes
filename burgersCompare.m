clear all
close all


nx = 1000;
dx = 2.*pi/(nx);
x = linspace(0,2*pi - dx,nx);

F = buildF(nx,dx);
G = inv(F);
u = initialCondition(x,5);

%% Get Fourier Coefficients and wave numbers
uhat = F*u;
k = linspace(-nx/2,nx/2 - 1,nx)';

t = 0;
dt = 0.01*dx;
et = 10*pi;
nu = 0.05;
iter = 0;
rk4const = [1./4,1./3,1./2,1.]
u2 = u
uhat2 =  uhat
while (t < et)
  iter = iter + 1;
  uhat0 = uhat;
  uhat02 = uhat2;
  for L=1:4
    RHS1 = convolutionRHS(uhat,k,nu);
    RHS2 = pseudoSpectralRHS(uhat2,k,nu,F,G); 
    %ur = G*uhat;
    %u2 =ur.*ur;
    %ghat = F*u2;
    uhat = uhat0 + dt*rk4const(L)*(RHS1);
    uhat2 = uhat02 + dt*rk4const(L)*RHS2;
  end
  t = t + dt
  if (mod(iter,10) == 0)
    subplot(1,2,1)
    plot(x,real(G*uhat))
    hold on
    plot(x,real(G*uhat2))
    hold off
    xlabel('x')
    ylabel('u')
    subplot(1,2,2)   
    semilogy(k,abs(uhat))
    %loglog(k(end/2+2:end),abs(uhat(end/2+2:end)))
    hold on
    semilogy(k,abs(uhat2))
    %loglog(k(end/2+2:end),abs(uhat2(end/2+2:end)))
    xlabel('k')
    ylabel('uhat')
    ylim([1e-10,4])
    %xlim([1,max(k)])
    hold off
    drawnow
  end
end
