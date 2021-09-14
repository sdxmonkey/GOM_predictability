% test the MC_maxchance,MC_minchance and MC_combinechance test
clear
ptest = [85,90,95,99];
Nrepeat = 10;
Nsample = 155;
Ntest = 100000;
Nevent = 5;
Ncombine = 3;

PN_max = nan([length(ptest),Nrepeat]);
PN_min = nan([length(ptest),Nrepeat]);
PN_comb = nan([2,length(ptest),Nrepeat]);
for irepeat = 1:Nrepeat
    irepeat
    tic
    % test MC_maxchance
    [nn,pn] = MC_maxchance(Nsample,Ntest,Nevent,ptest);
    if irepeat <= 6
        figure(1)
        subplot(2,3,irepeat)
        hist(nn)
        hold on
        for ii = 1:length(ptest)
            plot([pn(ii),pn(ii)],[0 6e4],'r')
        end
        hold off
        title(['MC_max ',num2str(irepeat)])
    end
    PN1(:,irepeat) = pn;
    % test MC_minchance
    [nn,pn] = MC_minchance(Nsample,Ntest,Nevent,ptest);
    if irepeat <= 6
        figure(2)
        subplot(2,3,irepeat)
        hist(nn)
        hold on
        for ii = 1:length(ptest)
            plot([pn(ii),pn(ii)],[0 6e4],'r')
        end
        hold off
        title(['MC_min ',num2str(irepeat)])
    end
    PN2(:,irepeat) = pn;
    % test MC_combinechance
    [nn,pn] = MC_combinechance(Nsample,Ntest,Nevent,Ncombine,ptest);
    if irepeat <= 6
        figure(3)
        subplot(2,3,irepeat)
        hist(nn)
        hold on
        for ii = 1:length(ptest)
            plot([pn(1,ii),pn(1,ii)],[0 6e4],'r')
            plot([pn(2,ii),pn(2,ii)],[0 6e4],'r')
        end
        hold off
        title(['MC_combine ',num2str(irepeat)])
    end
    PN3(:,:,irepeat) = pn;
    toc
end

figure(4);clf
for iplot = 1:4
    subplot(2,2,iplot)
    hist(PN1(iplot,:))
    title(['MC_max: ',num2str(ptest(iplot)),'%'])
end

figure(5);clf
for iplot = 1:4
    subplot(2,2,iplot)
    hist(PN2(iplot,:))
    title(['MC_min: ',num2str(ptest(iplot)),'%'])
end

figure(6);clf
p_bi = icdf('Binomial',ptest/100,Nsample,1/Nevent*Ncombine);
for iplot = 1:4
    subplot(2,2,iplot)
    hist(squeeze(PN3(1,iplot,:)))
    title(['MC_combine high: ',num2str(ptest(iplot)),'%'])
    hold on
    plot([p_bi(iplot) p_bi(iplot)],[0 Nrepeat],'r')
    hold off
end

figure(7);clf
p_bi = icdf('Binomial',1-ptest/100,Nsample,1/Nevent*Ncombine);
for iplot = 1:4
    subplot(2,2,iplot)
    hist(squeeze(PN3(2,iplot,:)))
    title(['MC_combine high: ',num2str(ptest(iplot)),'%'])
    hold on
    plot([p_bi(iplot) p_bi(iplot)],[0 Nrepeat],'r')
end
