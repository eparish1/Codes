
function [Flux] = convolutionRHS(uhat,k,nu,F,G)
  ur = G*uhat; % perform ifft
  u2 =ur.*ur; %square real component
  ghat = F*u2; %take FFT
  Flux = -0.5.*1j.*k.*ghat - nu.*k.^2.*uhat;

