function [J] = objectiveFunction(u,beta,d,beta0)
  J = (d - u).'*1e2*(d - u) + (beta - beta0).'*0*(beta - beta0);
end
