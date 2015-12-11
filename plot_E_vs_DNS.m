function [] = plot_E_vs_DNS(k,kc,uhatsave,DNS,GfDNS,tsave,t,string,print)    
    close all
    [t1,DNSind] = min(abs(DNS.tsave - t));
    [t1,LESind] = min(abs(tsave - t));
    uhat = uhatsave(:,LESind);
    loglog(abs(k(end/2+2:end)),0.5.*uhat(end/2+2:end).*conj(uhat(end/2+2:end)),'LineWidth',2)
    hold on
    loglog(abs(DNS.k),GfDNS.*DNS.uhatsave(:,DNSind).*conj(GfDNS.*DNS.uhatsave(:,DNSind)),'color','red','LineWidth',2)
    plot(ones(2,1)*kc,linspace(1e-10,1,2),'color','black','LineWidth',3)
    %loglog(abs(k),abs(k).^(-5./3.)*3.)
    %semilogy(k,abs(c).*exp(-k.^2*(50.*dx).^2./24.))
    xlabel('$k$','Interpreter','LaTex','FontSize',20)
    ylabel('$E(k)$','Interpreter','LaTex','FontSize',20)
    legend('LES','Filtered DNS','Filter Cutoff','Location','southeast')
    ylim([1e-10,40])
    xlim([1,5.*max(k)])
    if (print == 1)
      saveas(gcf,string)
    end
