function [tauhat] = computeSGS(Gf,uhat)
  tauhat = 0.5*(Gf.*myconv(uhat,uhat,nx) - myconv(Gf.*uhat,Gf.*uhat,nx));
end
