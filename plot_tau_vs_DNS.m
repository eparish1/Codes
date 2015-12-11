function [] = plot_tau_vs_DNS(k,kc,tauhatsave,DNS,GfDNS,tsave,t,string,print)
    close all
    [t1,DNSind] = min(abs(DNS.tsave - t));
    [t1,LESind] = min(abs(tsave - t));
    tauhat = tauhatsave(:,LESind);
    loglog(k,abs(tauhat),'linewidth',2)
    hold on
    loglog(DNS.k,abs(computeSGS(GfDNS,DNS.uhatsave(:,DNSind))),'color','red','linewidth',2)
    plot(ones(2,1)*kc,linspace(1e-10,1,2),'color','black','LineWidth',3)
    hold off
    %hold on
    %plot(DNS.x,real(ifft(fftshift(GfDNS.*DNS.uhatsave(:,DNSind)*1000))),'color','red')
    xlabel('$k$','Interpreter','LaTex','FontSize',20)
    ylabel('$|\tau_{sgs}|$','Interpreter','LaTex','FontSize',20)
    legend('LES','Filtered DNS','Filter Cutoff','Location','southeast')
    xlim([1,5*max(k)])
    if (print == 1)
      saveas(gcf,string)
    end
end
