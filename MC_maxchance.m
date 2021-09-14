function [nn,pn] = MC_maxchance(Nsample,Ntest,Nevent,ptest)
%MC_MAXCHANCE test the significance level that a target event has the
%largest chance
%output:
%       nn: The number of occurrences for the target event in each test [1xNtest]
%       pn: The number of occurrences for the target event at selected
%       significance level [1xNp]
%input:
%       Nsampe: sample size in each test.
%       Ntest: the number of times to repeat the test
%       Nevent：the total types of events
%       ptest: the siginificance level to test [1xNp]
%by sdxmonkey on Aug 8, 2021
nn = nan(1,Ntest);
P = nan(1,Nevent);
for itest = 1:Ntest
    % randomly generate the chance for each event
%     P(1) = rand(1);
%     Prest = 1-P(1);
%     for ievent = 2:Nevent-1
%         P(ievent) = Prest*rand(1);
%         Prest = Prest - P(ievent);
%     end
%     P(end) = Prest;
    tmp = sort([0,rand(1,4),1]);
    P = diff(tmp);
    P = sort(P);
    Ptarget = P(randi(4)); % target event will never be the event with highest chance
    
    xtest = rand([1,Nsample]);
    II = find(xtest <= Ptarget);
    nn(itest) = length(II);
end
pn = prctile(nn,ptest);
end

