function [u] = myifft(uhat,nx)
  u = (ifft(fftshift(uhat)))*nx;
end
