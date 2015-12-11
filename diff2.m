function [uxx] = diff2(u,dx)
  uxx = zeros(size(u));
  uxx(2:end-1) = (u(3:end) - 2.*u(2:end-1) + u(1:end-2))./dx.^2;
  uxx(end) = (u(1) - 2.*u(end) + u(end-1))./dx^2;
  uxx(1) =  (u(2) - 2.*u(1) + u(end))./dx^2;
end
