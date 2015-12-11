function [Q] =   advanceSolutionUnsteadyAdjoint(Q,beta)
  global turbmodel Cs Delta tausgs dt
  rk4const = [1./4,1./3,1./2,1.];
  Q0 = Q;
  for L = 1:4
    [RHS] = burgersInferRHS(Q,beta);
    Q = Q0 + dt*rk4const(L)*RHS;
  end
end
  

    
