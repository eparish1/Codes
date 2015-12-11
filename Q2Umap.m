function  Q2Umap()
  global turbmodel u uhat tausgs tauhatsgs Q solution_domain w0hat w1hat w2hat
  if (solution_domain == 1 )
    if (turbmodel <= 1 || turbmodel == 5)
      u = Q;
    end
    if (turbmodel >= 2 && turbmodel <=4)
      u = Q(1:2:end);
      tausgs = Q(2:2:end);
    end
  end
  if (solution_domain == 2)
    if (turbmodel <= 10)
      uhat = Q;
    end
    
    if (turbmodel == 11 || turbmodel == 12)
      uhat  = Q(1:3:end);
      w0hat = Q(2:3:end);
      w1hat = Q(3:3:end);
    end
    if (turbmodel == 13 || turbmodel == 14)
      uhat  = Q(1:4:end);
      w0hat = Q(2:4:end);
      w1hat = Q(3:4:end);
      w2hat = Q(4:4:end);
    end

  end
end
