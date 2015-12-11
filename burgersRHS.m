function [RHS] = burgersRHS(Q)
  global turbmodel k Cs Delta tausgs tauhatsgs nx nu dx u convective_flux solution_domain ...
         uhat M t x Gf w0hat w1hat w2hat
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
    if (turbmodel == 5)
      %kf = linspace(2*min(k),2*max(k+1)-1,2*nx);
      %GfD = heaviside(pi/Delta - abs(kf) )';
      %uhat = myifft(Q,nx);
      %Lhat = -0.5*cconv([uhat;real(uhat(1)) - 1j*imag(uhat(1))],[uhat;real(uhat(1)) - 1j*imag(uhat(1))],2*nx+1);
      %Lhat(nx/2+1:3*nx/2+1) = 0.;
      %dxf = 2.*pi/(nx*2);
      %xf = linspace(0,2.*pi - dxf,2*nx+1);
      %L = interp1(xf,real(myifft(Lhat,2.*nx)),x)';
      L = 0.5*smooth(u.*u,'lowess') - 0.5*smooth(u,'lowess').*smooth(u,'lowess');
      tausgs = -t.*u.*diff1(L,dx)*0.1;
      tauhatsgs = myfft(tausgs,nx);
      RHS = -f + nu.*uxx -diff1(tausgs,dx);
    end
  end
  if (solution_domain == 2)
    if (turbmodel == 0)
      c = myconv(Q,Q,length(Q));
      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*Q;
      tauhatsgs = 0.5*(Gf.*myconv(Q,Q,nx) - myconv(Gf.*Q,Gf.*Q,nx));
      tausgs = myifft(tauhatsgs,nx);
    end
    if (turbmodel == 1)
      utmp = myifft(uhat,nx);
      [nut,tau_eq] = smagorinsky(utmp,dx);
      nut_x = diff1(nut,dx);
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
       tauhatsgs = .35*myconv(uhat,tmp,length(uhat));
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
    if (turbmodel == 7)
       uhat = Q;
       DeltaD = Delta*2;
       GfD = heaviside(pi/DeltaD - abs(k));
       uhatFilt = GfD.*uhat;
       uReal = myifft(uhat,nx);
       uReal_x = diff1(uReal,dx);
       uRealFilt = myifft(uhatFilt,nx);
       uRealFilt_x =diff1(uRealFilt,dx);
       L = 0.5*GfD.*myconv(uhat,uhat,nx) - 0.5*myconv(GfD.*uhat,GfD.*uhat,nx);
       val = myfft(Delta^2*abs(uReal_x).*uReal_x,nx);
       val2 = myfft(DeltaD^2*abs(uRealFilt_x).*uRealFilt_x,nx);
       M = GfD.*val - val2;
       Cs2 = L./M;
       nut = Cs2.*Delta.^2.*abs(uReal_x);
       tausgs = -nut.*uReal_x;
       tauhatsgs = myfft(tausgs,nx);
       c = myconv(uhat,uhat,length(uhat));
       RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhatsgs;
    end
    if (turbmodel == 8)
      DeltaD = Delta*2;
      GfD = heaviside(pi/DeltaD - abs(k));
      L = 0.5*GfD.*myconv(Q,Q,nx) - 0.5*myconv(GfD.*Q,GfD.*Q,nx);
      tauhatsgs = 0.25*L;
      tausgs = myifft(tauhatsgs,nx);
      c = myconv(Q,Q,length(Q));
      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*Q - 1j.*k.*tauhatsgs;
    end

    if (turbmodel == 9)
      uhat2 = zeros(2*nx,1);
      uhat2(nx/2+1:nx/2+nx,1) = Q;
      uhat2(nx/2+nx+1) = real(Q(1)) - 1j*imag(Q(1));
      kf = linspace(-nx,nx-1,2*nx)';
      p = zeros(2*nx,1);
      q = kf - p;
      p(nx/2+1:nx/2+nx+1) = kf(nx/2+1:nx/2+nx+1);
      res = zeros(2*nx,1); res(nx/2+1:nx/2+nx+1) = 1;
      unres = ones(2*nx,1) - res;
      nxf = length(uhat2);
      PLu = -1j.*kf./2.*myconv(uhat2,uhat2,nxf) - nu.*kf.^2.*uhat2;
      PLu_p = res.*PLu;
      PLu_q = unres.*PLu;

       tmp = -t.*1j.*M./2.*myconv(uhat2,uhat2,length(uhat2));
       tauhatsgs = 0.35*myconv(uhat,tmp,length(uhat));

      PLPLu = -1j.*kf./2.*myconv(2.*uhat2,-1j.*p./2.*myconv(uhat2,uhat2,nxf) - nu.*p.^2.*uhat2 ,nxf) ...
              - p.^2.*nu.*(-1j.*kf./2.*myconv(uhat2,uhat2,nxf) - kf.^2.*nu.*uhat2);
      PLPLu_p = res.*PLPLu;
      PLPLu_q = unres.*PLPLu;
%      tau1 = t*myconv(uhat2,PLu_q,nxf);
%      tau2 = -t.^2/2.*(myconv(uhat2,PLPLu_q,nxf) + myconv(PLu_p,PLu_q,nxf));
%      tauhatsgs =  (tau1 + 1.0*tau2);
%            
      tmod = t.*2.*(-1j.*kf./2.*myconv(uhat2,PLu_q,nxf) - 0.5.*q.^2.*nu.*PLu_q);
      t2mod = -t.^2./2.*2.*(-1j.*kf./2.*myconv(uhat2,PLPLu_q,nxf) -0.5.*q.^2.*nu.*PLPLu_q ...
             +  -1j.*kf./2.*myconv(PLu_p,PLu_q,nxf) -0.5.*p.^2.*nu.*PLu_q);
      tterm = tmod(nx/2+1:nx/2+nx)+ t2mod(nx/2+1:nx/2+nx);
      tauhatsgs = tterm./(-1j.*k + 1e-30);
%      t2mod = -1j.*kf.*(myconv(uhat2,PLPLu_q,nxf)  + myconv(PLu_p,PLu_q,nxf));
%      tauhatsgs = t*myconv(uhat2,PLu_q,nxf) - t.^2/2.*(myconv(uhat2,PLPLu_q,nxf)  + myconv(PLu_p,PLu_q,nxf));
%      tauhatsgs = tauhatsgs(nx/2+1:nx/2+nx);
      tausgs = myifft(tauhatsgs,nx);

      c = myconv(uhat,uhat,length(uhat));
%      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*Q -1j.*k.*tauhatsgs;
      RHS =  -0.5.*1j.*k.*c - nu.*k.^2.*Q + tterm;
    end

    if (turbmodel == 10)
      L = computeLeonardStress(Gf,Q);
      tauhatsgs = 0.25*L;
      tausgs = myifft(tauhatsgs,nx);
      c = myconv(Q,Q,length(Q));
      RHS = -0.5.*1j.*k.*c - nu.*k.^2.*Q -1j.*k.*tauhatsgs;
    end
    if (turbmodel == 11)
      Q2Umap();
      uhat2 = zeros(2*nx,1);
      kf = linspace(-nx,nx-1,2*nx)';
      uhat2(nx/2+1:nx/2+nx,1) = uhat;
      uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
      nxf = length(uhat2);
      PLu = -1j.*kf./2.*myconv(uhat2,uhat2,nxf) - nu.*kf.^2.*uhat2;
      PLu_q = PLu; PLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLu_p = PLu - PLu_q;

      PLQLu = 2.*-1j.*kf./2.*myconv(uhat2,PLu_q,nxf) - kf.^2.*nu.*PLu_q; 
      PLQLu_q = PLQLu; PLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLu_p = PLQLu - PLQLu_q;

      PLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLu_q,uhat2) + 2.*-1j.*kf./2.*myconv(PLu_q,PLu_p) + ...
                  2.*-1j.*kf./2.*myconv(PLu_q,PLu_q) - nu.*kf.^2.*PLQLu_q;

      RHS1 = PLu(nx/2+1:nx/2+nx) + 0.1*w0hat;
      RHS2 = PLQLu(nx/2+1:nx/2+nx) + 0.0*w1hat;
      RHS3 = PLQLQLu(nx/2+1:nx/2+nx) - 100.*k.^2*nu.*w1hat;
      tauhatsgs = 0.02*w0hat./(-1j.*k + 1e-30);
      tausgs = myifft(tauhatsgs,nx);
      RHS = zeros(nx*3,1);
      RHS(1:3:end) = RHS1;
      RHS(2:3:end) = RHS2;
      RHS(3:3:end) = RHS3; 
    end
    if (turbmodel == 12)
      dt0 = 0.1;
      dt1 = 0.01;
      Q2Umap();
      uhat2 = zeros(2*nx,1);
      kf = linspace(-nx,nx-1,2*nx)';
      uhat2(nx/2+1:nx/2+nx,1) = uhat;
      uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
      nxf = length(uhat2);
      PLu = -1j.*kf./2.*myconv(uhat2,uhat2,nxf) - nu.*kf.^2.*uhat2;
      PLu_q = PLu; PLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLu_p = PLu - PLu_q;

      PLQLu = 2.*-1j.*kf./2.*myconv(uhat2,PLu_q,nxf) - kf.^2.*nu.*PLu_q; 
      PLQLu_q = PLQLu; PLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLu_p = PLQLu - PLQLu_q;

      PLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLu_q,uhat2) + 2.*-1j.*kf./2.*myconv(PLu_q,PLu_p) + ...
                  2.*-1j.*kf./2.*myconv(PLu_q,PLu_q) - nu.*kf.^2.*PLQLu_q;

      RHS1 = PLu(nx/2+1:nx/2+nx) + w0hat;
      RHS2 = -2./dt0.*w0hat + 2.*PLQLu(nx/2+1:nx/2+nx) + w1hat;
      RHS3 = -2./dt1.*w1hat + 2.*PLQLQLu(nx/2+1:nx/2+nx);
      tauhatsgs = w0hat./(-1j.*k + 1e-30);
      tausgs = myifft(tauhatsgs,nx);
      RHS = zeros(nx*3,1);
      RHS(1:3:end) = RHS1;
      RHS(2:3:end) = RHS2;
      RHS(3:3:end) = RHS3; 
    end


    if (turbmodel == 13)
      Q2Umap();
      uhat2 = zeros(2*nx,1);
      kf = linspace(-nx,nx-1,2*nx)';
      uhat2(nx/2+1:nx/2+nx,1) = uhat;
      uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
      nxf = length(uhat2);
      PLu = -1j.*kf./2.*myconv(uhat2,uhat2,nxf) - nu.*kf.^2.*uhat2;
      PLu_q = PLu; PLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLu_p = PLu - PLu_q;


      PLQLu = 2.*-1j.*kf./2.*myconv(uhat2,PLu_q,nxf) - kf.^2.*nu.*PLu_q; 
      PLQLu_q = PLQLu; PLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLu_p = PLQLu - PLQLu_q;

      PLPLu  = 2.*-1j.*kf./2.*myconv(uhat2,PLu_p) - kf.^2.*nu.*PLu_p;
      PLPLu_q = PLPLu; PLPLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLPLu_p = PLPLu - PLPLu_q;


      PLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLu_q,uhat2) + 2.*-1j.*kf./2.*myconv(PLu_q,PLu_p) + ...
                  2.*-1j.*kf./2.*myconv(PLu_q,PLu_q) - nu.*kf.^2.*PLQLu_q;

      PLQLQLu_q = PLQLQLu; PLQLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLQLu_p = PLQLQLu - PLQLQLu_q;

      PLQLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLQLu_q,uhat2) + 4.*-1j.*kf./2.*myconv(PLu_p,PLQLu_q) + ...
                  4.*-1j.*kf./2.*myconv(PLQLu_p,PLu_q) + 2.*-1j.*kf./2.*myconv(PLPLu_p,PLu_q) + ...
                  6.*-1j.*kf./2.*myconv(PLQLu_q,PLu_q) + 2.*-1j.*kf./2.*myconv(PLPLu_q,PLu_q) + ...
                  - nu.*kf.^2.*PLQLQLu_q;

   
      RHS1 = PLu(nx/2+1:nx/2+nx) + w0hat;
      RHS2 = PLQLu(nx/2+1:nx/2+nx) + w1hat;
      RHS3 = PLQLQLu(nx/2+1:nx/2+nx) + w2hat;
      dt2 = 1;
      RHS4 = -2./dt2.*w2hat + 2.*PLQLQLQLu(nx/2+1:nx/2+nx);
      tauhatsgs = w0hat./(-1j.*k + 1e-30);
      tausgs = myifft(tauhatsgs,nx);
      RHS = zeros(nx*4,1);
      RHS(1:4:end) = RHS1;
      RHS(2:4:end) = RHS2;
      RHS(3:4:end) = RHS3; 
      RHS(4:4:end) = RHS4; 
    end



    if (turbmodel == 14)
      dt0 = 0.13;
      dt1 = 0.07;
      dt2 = 0.05;
      Q2Umap();
      uhat2 = zeros(2*nx,1);
      kf = linspace(-nx,nx-1,2*nx)';
      uhat2(nx/2+1:nx/2+nx,1) = uhat;
      uhat2(nx/2+nx+1) = real(uhat(1)) - 1j*imag(uhat(1));
      nxf = length(uhat2);
      PLu = -1j.*kf./2.*myconv(uhat2,uhat2,nxf) - nu.*kf.^2.*uhat2;
      PLu_q = PLu; PLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLu_p = PLu - PLu_q;


      PLQLu = 2.*-1j.*kf./2.*myconv(uhat2,PLu_q,nxf) - kf.^2.*nu.*PLu_q; 
      PLQLu_q = PLQLu; PLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLu_p = PLQLu - PLQLu_q;

      PLPLu  = 2.*-1j.*kf./2.*myconv(uhat2,PLu_p) - kf.^2.*nu.*PLu_p;
      PLPLu_q = PLPLu; PLPLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLPLu_p = PLPLu - PLPLu_q;


      PLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLu_q,uhat2) + 2.*-1j.*kf./2.*myconv(PLu_q,PLu_p) + ...
                  2.*-1j.*kf./2.*myconv(PLu_q,PLu_q) - nu.*kf.^2.*PLQLu_q;

      PLQLQLu_q = PLQLQLu; PLQLQLu_q(nx/2+1:nx/2+nx+1) = 0.;
      PLQLQLu_p = PLQLQLu - PLQLQLu_q;

      PLQLQLQLu = 2.*-1j.*kf./2.*myconv(PLQLQLu_q,uhat2) + 4.*-1j.*kf./2.*myconv(PLu_p,PLQLu_q) + ...
                  4.*-1j.*kf./2.*myconv(PLQLu_p,PLu_q) + 2.*-1j.*kf./2.*myconv(PLPLu_p,PLu_q) + ...
                  6.*-1j.*kf./2.*myconv(PLQLu_q,PLu_q) + 2.*-1j.*kf./2.*myconv(PLPLu_q,PLu_q) + ...
                  - nu.*kf.^2.*PLQLQLu_q;

   
      RHS1 = PLu(nx/2+1:nx/2+nx) + w0hat;
      RHS2 = -2./dt0.*w0hat + 2.*PLQLu(nx/2+1:nx/2+nx) + w1hat;
      RHS3 = -2./dt1.*w1hat + 2.*PLQLQLu(nx/2+1:nx/2+nx) + w2hat;
      RHS4 = -2./dt2.*w2hat + 2.*PLQLQLQLu(nx/2+1:nx/2+nx);
      tauhatsgs = w0hat./(-1j.*k + 1e-30);
      tausgs = myifft(tauhatsgs,nx);
      RHS = zeros(nx*4,1);
      RHS(1:4:end) = RHS1;
      RHS(2:4:end) = RHS2;
      RHS(3:4:end) = RHS3; 
      RHS(4:4:end) = RHS4; 
    end

end
