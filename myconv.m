function [uwhatconv] = myconv(uhat,what,nx)
   uhat = [uhat;real(uhat(1,1))-1j*imag(uhat(1,1))];
   what = [what;real(what(1,1))-1j*imag(what(1,1))];
   uwhatconv = conv(uhat,what,'same');
   uwhatconv = uwhatconv(1:end-1);
   %uwhatconv(end) = real(uwhatconv(2)) - 1j.*imag(uwhatconv(2));
   %uwhatconv = fftshift(cconv(uhat,what,nx));
