clear all
close all


nx = 500;
L = 2.*pi;
dx = L/(nx);
x = linspace(0,2.*pi - dx,nx);

k = linspace(-nx/2,nx/2 - 1,nx)';
u = sin(x)' + 4;
%% Get Fourier Coefficients and wave numbers
%uhat = F*u;
uhat = myfft(u,nx);
h = figure;

u2 = u;
t = 0;
et = 5;
nu = 5e-2;
%dt = nu*dx^2*20;
dt = 2e-3;
iter = 0;
rk4const = [1./4,1./3,1./2,1.];
iter2 = 1;
%%% Filter Info
Delta = 50.*dx;

while (t < et)
  
  iter = iter + 1;
  uhat0 = uhat;
  u20 = u2;
  K0 = 0.5.*u2.*u2;
  for L=1:1
    [RHS1,c] = convolutionRHS(uhat,k,nu);
    Jac = admDiffFD(@convolutionRHS,1,uhat,k,nu,admOptions('i',1,'d',1));
    A = eye(nx)./dt - Jac;
    uhat = A\RHS1 + uhat;
    [RHS2] = upwindRHS(u2,nu,dx);
    Jac2 = admDiffFD(@upwindRHS,1,u2,nu,dx,admOptions('i',1,'d',1));
    A2 = eye(nx)./dt - Jac2;
    u2 = A2\RHS2 + u2;
    %uhat = uhat0 + dt*(RHS1);
    %u2 = u20 + dt*RHS2;
  end
  t = t + dt;
  K = 0.5*u2.*u2;
  saveval(iter) = norm(abs(1./dt*(uhat - uhat0) - RHS1));
  saveval2(iter) = norm(1./dt*(K - K0) + 2./3.*diff1(u20.*K0,dx) - u*nu.*diff2(u,dx) );
  if (mod(iter,10) == 0)
    iter2 = iter2 + 1;
    subplot(1,2,1)
    plot(x,real(myifft(uhat,nx)))
    hold on
    plot(x,u2,'o','color','red')
    drawnow
    hold off
    t
  end
end



