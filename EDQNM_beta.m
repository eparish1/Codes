function [tau_hat] = EDQNM_beta(k,kc,uhat,beta)
  [val,index] = min(abs(k - kc));
  nut = real(beta).*(uhat(index).*conj(uhat(index))/kc).^real(0.35);
  tau_hat = -1j.*k.*nut.*uhat;
end

