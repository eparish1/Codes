%function [u] = initialCondition(x,N,k)
%  u = zeros(size(x));
%  u(:) = 4;
%  for i = 1:N
%    u = u +  sin(x*i + i^2);
%  end
%  u = u'; 

function [u] = initialCondition(x,kc,k)
  u = zeros(size(x));
  E = zeros(size(x));
  rng(1); %seed random number generator
  beta = rand(1,length(x))*2*pi - pi;
  for i = length(x)/2+2:length(x)
    if  ( k(i) >=1) & ( k(i) <=5)
      E(i) = 5^(-5./3.);
    else
      if k(i) < kc
      E(i) = abs(k(i))^(-5./3.);
      end
    end
    u = u + sqrt(2.*E(i)).*sin(abs(k(i)).*x + beta(i));
  end
  u = u';
