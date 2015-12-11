function  U2Qmap()
  global turbmodel u uhat tausgs Q nx solution_domain tauhatsgs w0hat w1hat w2hat 
  if (solution_domain == 1)
    if (turbmodel <= 1 || turbmodel == 5)
      Q = u;
    end
    if (turbmodel >= 2 && turbmodel <= 4)
      Q = zeros(2*nx,1);
      Q(1:2:end) = u;
      Q(2:2:end) = tausgs;
    end
  end

  if (solution_domain == 2)
    if (turbmodel <= 10)
      Q = uhat;
    end
    if (turbmodel == 11 || turbmodel == 12)
      Q = zeros(3*nx,1);
      Q(1:3:end) = uhat;
      Q(2:3:end) = w0hat;
      Q(3:3:end) = w1hat;
    end  
    if (turbmodel == 13 || turbmodel == 14)
      Q = zeros(4*nx,1);
      Q(1:4:end) = uhat;
      Q(2:4:end) = w0hat;
      Q(3:4:end) = w1hat;
      Q(4:4:end) = w2hat;
    end  


  end
