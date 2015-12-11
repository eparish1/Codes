function  advanceSolutionInfer()
  global turbmodel Cs Delta tausgs dt nx Q beta dQdBeta
  rk4const = [1./4,1./3,1./2,1.];
  Q0 = Q;
  dQdBeta0 = dQdBeta;
  for L = 1:1
    [RHS] = burgersInferRHS(Q,beta);
    dhoRdhoQ = admDiffFD(@burgersInferRHS,1,Q,beta,admOptions('i',1,'d',1));
    dhoRdhoBeta =  admDiffFD(@burgersInferRHS,1,Q,beta,admOptions('i',2,'d',1));
    dQdBeta = dQdBeta0 + dt.*((dhoRdhoQ*dQdBeta) - dhoRdhoBeta);
    Q = (eye(nx)./dt - dhoRdhoQ)\RHS + Q0;
  end
end
  

    
