function [F] = kurgTadmorFlux(u,nu,dx)
  phi = fluxLim(u);
  [uLp,uLm,uRp,uRm] = MUSCL_LR(u,phi);
  Fp = zeros(size(u));
  Fm = zeros(size(u));
  a1 = zeros(size(u));
  uabs = abs(u);
  a1(2:end) = max(uabs(2:end),uabs(1:end-1));
  a1(1) = max(uabs(1),uabs(end));
  a2 = zeros(size(u));
  a2(1:end-1) = max(uabs(2:end),uabs(1:end-1));
  a2(end) = max(uabs(1),uabs(end));
  Fm = 0.5.*( (burgerFlux(uRm) + burgerFlux(uLm)) - a1.*(uRm - uLm));
  Fp = 0.5.*( (burgerFlux(uRp) + burgerFlux(uLp)) - a2.*(uRp - uLp));
  F = 1./dx.*(Fp - Fm);
end


function [F] = burgerFlux(u)
  F = 0.5.*u.^2;
end

function [uLp,uLm,uRp,uRm] =  MUSCL_LR(u,phi)
  uLp = zeros(size(u));
  uRp = zeros(size(u));
  uRm = zeros(size(u));
  uLm = zeros(size(u));
  uLp(1:end-1) = u(1:end-1) + 0.5.*phi(1:end-1).*(u(2:end) - u(1:end-1));
  uLp(end) = u(end) + 0.5*phi(end)*(u(1) - u(end));
  uLm(2:end) = u(1:end-1) + 0.5.*phi(1:end-1).*(u(2:end) - u(1:end-1));
  uLm(1) = u(end) + 0.5.*phi(end).*(u(1) - u(end));

  uRp(1:end-2) = u(2:end-1) - 0.5.*phi(2:end-1).*(u(3:end) - u(2:end-1));
  uRp(end-1) = u(end) - 0.5.*phi(end).*(u(1) - u(end));
  uRp(end) = u(1) - 0.5.*phi(1).*(u(2) - u(1));
  uRm(1:end-1) = u(1:end-1) - 0.5.*phi(1:end-1).*(u(2:end) - u(1:end-1));
  uRm(end) = u(end) - 0.5.*phi(end).*(u(1) - u(end));
end


function [phi] = fluxLim(u)
  r = zeros(size(u));
  r(2:end-1) = (u(2:end-1) - u(1:end-2))./(u(3:end) - u(2:end-1) + 1.e-20);
  r(1) = (u(1) - u(end))/(u(2) - u(1) + 1.e-20);
  r(end) = (u(end) - u(end-1))/(u(1) - u(end) + 1.e-20);
  phi = (r.^2 + r)./(r.^2 + 1);
end

