function [tau] = initialCondition(uhat,k)
  c1 = conv(uhat,uhat,'same');
  Ghat = exp(-k);
  c2 = conv(Ghat.*uhat,Ghat.*uhat);
  tau = zeros(size(uhat));
  tau = 0.5+Ghat*c1 - 0.5*c2
