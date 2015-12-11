function [uhat] = myfft(u,nx)
  uhat = fftshift(fft(u))/nx;
end
