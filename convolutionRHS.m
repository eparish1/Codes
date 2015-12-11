
function [Flux,c] = convolutionRHS(uhat,k,nu)
  %c = 1./length(uhat)*conv(uhat,uhat,'same');
  %c = 1./length(uhat)*fftshift(cconv(uhat,uhat,length(uhat)));
  c = myconv(uhat,uhat,length(uhat));
  Flux = -0.5.*1j.*k.*c - nu.*k.^2.*uhat;

