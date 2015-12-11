function [nut,tau_eq] = smagorinksy(u,dx)
  global  Cs Delta
  nut = zeros(size(u));
  ux = diff1(u,dx);
  nut = (Cs*Delta).^2.*abs(ux);
  tau_eq = -nut.*ux;
end


