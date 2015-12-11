function [Q,Qsave] = solverBurgers(u,uhat,t,et,save_freq,live_plot)
  tsave = [t];
  iter = 1;
  iter2 = 1;
  while (t < et)
    iter = iter + 1;
    advanceSolution();
    Q2Umap();
    t = t + dt;
    if (mod(iter,save_freq) == 0)
      [val,DNSind] = min(abs(t - DNS.tsave));
      tsave = [tsave,t];
      usave = [usave,u];
      uhatsave = [uhatsave;uhat];
      iter2 = iter2 + 1;
      if (live_plot == 1)
        subplot(1,2,1)
        plot(x,real(u))
        hold on
        plot(DNS.x,DNS.ufiltsave(:,DNSind));
        drawnow
        hold off
        subplot(1,2,2)
        plot(x,real(tausgs))
        hold on
        %plot(DNS.x,real(DNS.tauSGSsave(:,DNSind)))
        plot(x,interp1(DNS.x,DNS.tauSGSsave(:,DNSind),x))
        hold off
      end
      t
    end
  end
end
