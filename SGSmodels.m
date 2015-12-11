function [RHStauSGS] = RHSSGS(tauSGS,uhat,k)
   T1 = 2./3.*1j.*k.*conv(uhat,tauSGS,'same') ;  
   T2 = -k.^2*nu.*tauSGS;
   T3 = conv(uhat,1j.*k.*tauSGS,'same');
   %phi = 2./6.*(1j.*k.*conv(uhat,Gf.*conv(uhat,uhat,'same'),'same') - ...
   %      1j.*k.*Gf.*conv(uhat,conv(uhat,uhat,'same'),'same')) + ...
   %      nu.*conv(1j.*k.*ufilt,1j.*k.*ufilt,'same') - ...
   %      nu.*Gf.*(conv(1j.*k.*uhat,1j.*k.*uhat,'same'));
   RHStauSGS = T1 - T2 - T3;
end 

function [tauSGS] = smag(uhat,k,dx,Delta)
  u = ifft(fftshift(uhat));
  u_x = zeros(size(u));
  u_x(2:end-1) = (u(3:end) - u(1:end-2))/(2.*dx);
  u_x(end) =  (u(1) - u(end-1))/(2.*dx);
  u_x(1) =  (u(2) - u(end))/(2.*dx);
  u_x_hat = fftshift(fft(abs(u_x)));
  tauSGS = -(C*Delta)^2*conv(u_x_hat,i.*k.*uhat,'same')
end

  
 