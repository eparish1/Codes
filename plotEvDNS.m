function [] = plot_tau_vs_DNS(k,tauhat,DNS,GfDNS)
    loglog(abs(k),abs(tauhat))
    hold on
    loglog((DNS.k),abs(computeSGS(GfDNS,DNS.uhatsave(:,DNSind))),'color','red')
    plot(ones(2,1)*40,linspace(1e-10,1,2),'color','black','LineWidth',3)
    hold off
    %hold on
    %plot(DNS.x,real(ifft(fftshift(GfDNS.*DNS.uhatsave(:,DNSind)*1000))),'color','red')
    xlabel('$k$','Interpreter','LaTex','FontSize',20)
    ylabel('$|\tau_{sgs}|$','Interpreter','LaTex','FontSize',20)
    legend('DNS','Filtered DNS','Filter Cutoff')
end
