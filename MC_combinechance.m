function [nn,pn] = MC_combinechance(Nsample,Ntest,Nevent,Ncombine,ptest)
%MC_MAXCHANCE test the significance level that a set of combined events is
%prefered or not prefered. The null hypothesis is that all events have
%equal chances
%output:
%       nn: The number of occurrences for the combined event in each test [1xNtest]
%       pn: The number of occurrences for the target combine event at selected
%       significance level [2xNp]: row 1 is the occurrence if the target
%       event is perfered, row 2 is the occurrence if the target
%       event is not perfered
%input:
%       Nsampe: sample size in each test.
%       Ntest: the number of times to repeat the test
%       Nevent：the total types of events
%       Ncombine: 
%       ptest: the siginificance level to test [1xNp]
%by sdxmonkey on Aug 25, 2021
nn = nan(1,Ntest);
pn = nan(2,length(ptest));
for itest = 1:Ntest
    XX = rand([Nevent,Nsample]);
    [M,I] = max(XX,[],1);
    Ievent = find(I<=Ncombine);
    nn(itest) = length(Ievent);
end
pn(1,:) = prctile(nn,ptest);
pn(2,:) = prctile(nn,100 - ptest);
end

