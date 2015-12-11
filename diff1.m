function [up] = diff1(u,dx)
  up = zeros(size(u));
  up(2:end) = 1./dx*(u(2:end) - u(1:end-1));
  up(1) = 1./dx*(u(1) - u(end));
end
