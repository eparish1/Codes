
function [F] = buildF(N,h)
  F = zeros(N,N);
  for m=1:N
    for n=1:N
      %F(m,n) = h/(2.*pi)*exp(-1j*(m - N/2.-1)*(n-1)*h);
      F(m,n) = 1/N*exp(-1j*(m - N/2.-1)*(n-1)*h);

    end
  end


