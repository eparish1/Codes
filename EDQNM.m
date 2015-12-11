function [tau_hat] = EDQNM(k,kc,uhat)
  [val,index] = min(abs(k - kc));
  nut = 0.28.*ones(size(k)).*sqrt(uhat(index).*conj(uhat(index))/kc);
  tau_hat = -1j.*k.*nut.*uhat;
end

