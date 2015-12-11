function [RHS,tauhat] = convolutionRHS_LES(uhat,tauhat,k,kc,nu)
  %c = 1./(length(uhat))*conv(uhat,uhat,'same');
  c = myconv(uhat,uhat,length(uhat));
  tauhat = EDQNM(k,kc,uhat);
  RHS = -0.5.*1j.*k.*c - nu.*k.^2.*uhat - 1j.*k.*tauhat;
end

