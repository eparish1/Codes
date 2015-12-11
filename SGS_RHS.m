
function [Flux] = convolutionRHS(uhat,k,nu)
  c = conv(uhat,tauhat,'same');
  g = conv(uhat,k*tauhat,'same');
  Flux = -0.5.*1j.*k.*c - nu.*k.^2.*tauhat + i*g;

