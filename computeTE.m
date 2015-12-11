function [TE] = computeTE(k,kc,uhat)
  [nx,nsaves] = size(uhat);
  TE = [];
  [dum,j_min] = min(abs(k - -kc));
  [dum,j_max] = min(abs(k - (kc-1)));
  for i=1:nsaves
    TE(i) = 0.5*sum(uhat(j_min:j_max,i).*conj(uhat(j_min:j_max,i)));
  end
  TE = TE;
end
