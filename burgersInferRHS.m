function [RHS] = burgersInferRHS(Q,beta)
  global turbmodel k Cs Delta tausgs tauhatsgs nx nu dx u convective_flux solution_domain  M t x Gf 
  if (solution_domain == 1)
    Q2Umap();
    uxx = zeros(size(u));
    uxx(2:end-1) = (u(3:end) - 2.*u(2:end-1) + u(1:end-2))./dx^2;
    uxx(end) = (u(1) - 2.*u(end) + u(end-1))./dx^2;
    uxx(1) =  (u(2) - 2.*u(1) + u(end))./dx^2;
    if (convective_flux == 1)
      %% Upwind Convective Flux
      fplus = 0.25*u.*(u + abs(u));
      fminus = 0.25*u.*(u - abs(u));
      f = zeros(size(u));
      f(2:end) = 1./dx*(fplus(2:end) - fplus(1:end-1));
      f(1) = 1./dx*(fplus(1) - fplus(end));
      f(1:end-1) = f(1:end-1) + 1./dx*(fminus(2:end) - fminus(1:end-1));
      f(end) = f(end) + 1./dx*(fminus(1) - fminus(end));
    end

    if (convective_flux == 2)
      f = kurgTadmorFlux(u,nu,dx);
    end
    if (convective_flux == 3)
      f = zeros(size(u));
      f(2:end-1) = (u(3:end).^2 - u(1:end-2).^2)/(2.*dx);
      f(end) = (u(1).^2 - u(end-1).^2)/(2.*dx);
      f(1) = (u(2).^2 - u(end).^2)/(2.*dx);
    end
    if (turbmodel == 0)
      RHS = -f + nu.*uxx;
      tausgs = zeros(size(u));
      tauhatsgs = zeros(size(u));
    end
    if (turbmodel == 1)
      [nut,tau_eq] = smagorinsky(u,dx);
      nut_x = diff1(nut,dx);
      u_x = diff1(u,dx);
      RHS = -f + nu.*uxx + nut_x.*u_x + nut.*uxx;
      tausgs = tau_eq;
      tauhatsgs = myfft(tausgs,nx);
    end
  
    if (turbmodel == 2)
      tscale = -0.1;
      u_x = diff1(u,dx);
      [nut,tau_eq] = smagorinsky(u,dx);
      tau_x = diff1(tausgs,dx);
      RHS1 = -f + nu.*uxx - tau_x;
      %%% Upwind tau equation
      ftplus = 0.5*(u + abs(u)).*tausgs;
      ftminus = 0.5*(u - abs(u)).*tausgs;
      ft = zeros(size(u));
      ft(2:end) = 1./dx*(ftplus(2:end) - ftplus(1:end-1));
      ft(1) = 1./dx*(ftplus(1) - ftplus(end));
      ft(1:end-1) = ft(1:end-1) + 1./dx*(ftminus(2:end) - ftminus(1:end-1));
      ft(end) = ft(end) + 1./dx*(ftminus(1) - ftminus(end));
      RHS2 =-ft + 1./tscale * (tausgs - tau_eq);
      RHS = zeros(2*nx,1);
      RHS(1:2:end) = RHS1;
      RHS(2:2:end) = RHS2;
    end

    if (turbmodel == 3)
      tscale = -0.2;
      u_x = diff1(u,dx);
      [nut,tau_eq] = smagorinsky(u,dx);
      tau_x = diff1(tausgs,dx);
      RHS1 = -f + nu.*uxx - tau_x;
      %%% Upwind tau equation
      uplus = 0.5*(u + abs(u));
      uminus = 0.5*(u - abs(u));
      ft = zeros(size(u));
      ft(2:end) = uplus(2:end).*1./dx.*(tausgs(2:end) - tausgs(1:end-1));
      ft(1) = uplus(1).*1./dx.*(tausgs(1) - tausgs(end));
      ft(1:end-1) = ft(1:end-1) + uminus(1:end-1).*1./dx.*(tausgs(2:end) - tausgs(1:end-1));
      ft(end) = ft(end) + uminus(end).*1./dx.*(tausgs(1) - tausgs(end));
      RHS2 =-ft + 1./tscale * (tausgs - tau_eq);
      RHS = zeros(2*nx,1);
      RHS(1:2:end) = RHS1;
      RHS(2:2:end) = RHS2;
    end

    if (turbmodel == 4)
      tau_x = diff1(tausgs,dx);
      [nut,tau_eq] = smagorinsky(u,dx);
      RHS1 = -f + nu.*uxx - tau_x;
      %%% Upwind tau equation
      ftplus = 0.5*(u + abs(u)).*tausgs;
      ftminus = 0.5*(u - abs(u)).*tausgs;
      ft = zeros(size(u));
      ft(2:end) = 1./dx*(ftplus(2:end) - ftplus(1:end-1));
      ft(1) = 1./dx*(ftplus(1) - ftplus(end));
      ft(1:end-1) = ft(1:end-1) + 1./dx*(ftminus(2:end) - ftminus(1:end-1));
      ft(end) = ft(end) + 1./dx*(ftminus(1) - ftminus(end));
      RHS = zeros(2*nx,1);
      RHS2 = -2./3.*ft + (nu + 15*nut).*diff2(tausgs,dx) + u.*diff1(tausgs,dx); 
      RHS(1:2:end) = RHS1;
      RHS(2:2:end) = RHS2;
    end
  end
  if (solution_domain == 2)
    if (turbmodel == 0)
      uhat = Q;
      c = myconv(uhat,uhat,length(uhat));
      c(nx/2+2:nx) = flipud(real(c(2:nx/2) )- 1j.*imag(c(2:nx/2)));
      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat;
      GfDNS = heaviside(pi/Delta - abs(k));
      tauhatsgs = 0.5*(GfDNS.*myconv(uhat,uhat,nx) - myconv(GfDNS.*uhat,GfDNS.*uhat,nx));
      tausgs = myifft(tauhatsgs,nx);
    end
    if (turbmodel == 1)
      uhat = Q;
      utmp = myifft(uhat,nx);
      [nut,tau_eq] = smagorinsky(utmp,dx);
      tausgs = tau_eq;
      tauhatsgs = myfft(tausgs,nx);
      c = myconv(uhat,uhat,length(uhat));
      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhatsgs;
    end
    if (turbmodel == 5)
       uhat = Q;
       c = myconv(uhat,uhat,length(uhat));
       uhat2 = zeros(2*nx,1);
       uhat2(nx/2+1:nx/2+nx,1) = uhat;
       uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
       tmp = -t.*1j.*M./2.*myconv(uhat2,uhat2,length(uhat2));
       tauhatsgs = real(beta).*myconv(uhat,tmp,length(uhat));
       RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhatsgs;
       tausgs = myifft(tauhatsgs,nx);
    end
    if (turbmodel == 6)
       kc = pi/Delta;
       c = myconv(uhat,uhat,length(uhat));
      [val,index] = min(abs(k - kc));
       nut = 0.47.*ones(size(k)).*(uhat(index).*conj(uhat(index))/kc).^0.37;
       tauhatsgs = -1j.*k.*nut.*uhat;
       %tauhatsgs = myconv(-nut,1j.*k.*uhat);
       RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhatsgs;
       tausgs = myifft(tauhatsgs,nx);
    end
  end
end
