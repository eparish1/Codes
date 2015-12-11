function [L] = computeLeonardStress(Gf,uhat)
  L = 0.5*Gf.*myconv(uhat,uhat,length(uhat)) - 0.5*myconv(Gf.*uhat,Gf.*uhat,length(uhat));
end


