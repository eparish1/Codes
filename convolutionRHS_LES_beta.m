function [RHS,tauhat] = convolutionRHS_LES_beta(uhat,k,kc,nu,beta)
  %c = 1./(length(uhat))*conv(uhat,uhat,'same');
  c = myconv(uhat,uhat,length(uhat));
  tauhat = EDQNM_beta(k,kc,uhat,beta);
  RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhat;
end

