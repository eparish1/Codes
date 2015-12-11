function  Q2Umap()
  global turbmodel u uhat tausgs Q solution_domain
  if (solution_domain == 1)
    if (turbmodel <= 1)
      u = Q;
    end
    if (turbmodel >= 2 && turbmodel <=4)
      u = Q(1:2:end);
      tausgs = Q(2:2:end);
    end
  end
  if (solution_domain == 2)
    if (turbmodel <= 1 || turbmodel == 5 || turbmodel == 6)
      uhat = Q;
    end
  end
end
