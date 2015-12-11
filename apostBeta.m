DNS = load('DNSdata/DNS.mat');
GfDNS = heaviside(20 - abs(DNS.k));%Sharp Spectral
for i=1:1001
  [val,index] = min(abs(DNS.k - kc));
  index = index ;
  ufilt(:,i) = GfDNS.*DNS.uhatsave(:,i);
  nut(:,i) = DNS.tauSGSsave(:,i)./(-1j.*DNS.k.*GfDNS.*ufilt(:,i));
  betaval1(:,i) = nut(end/2+2:end/2+21,i)./(  (ufilt(index,i).*conj(ufilt(index,i))./20).^real(0.35)  );
  ufilt(index,i)
end
