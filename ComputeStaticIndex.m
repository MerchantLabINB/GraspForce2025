function [StaticIndex,SigStaticIndex,PrepExeIndex] = ComputeStaticIndex(SigStaticClean,bootvalidationAccuracy,testvvalidationAccuracy)

Tbin = numel(SigStaticClean(:,1));

SigP = 0.05;

for bin = 1:Tbin  %%% over the trainin dimension compute the mean of the sig index in the testing dimansion
    Sboo = sort(bootvalidationAccuracy{bin});
    SigSboo = Sboo(950);
    sigvalidationAccuracy = 1:1:Tbin;
   % sigvalidationAccuracy = testvvalidationAccuracy(:,bin) > SigSboo;
    %StaticIndex(bin,1) = mean(SigStaticClean(sigvalidationAccuracy,bin));
    StaticIndex(bin,1) = mean(SigStaticClean(bin,sigvalidationAccuracy));
    [h,temtest] = ttest(SigStaticClean(sigvalidationAccuracy,bin))
    if temtest < SigP
        SigStaticIndex(bin,1) = 0.5;
    else
        SigStaticIndex(bin,1) = 0;
    end
    
end
SigStaticIndex(SigStaticIndex == 0) = NaN;

A = [1:61;1:61;1:61;1:61];
sz = size(A);
X = NaN(sz);
PrepExeIndex = X';

for bin = 1:Tbin%%% over the trainin dimension compute     
    PrepIndexPrep(bin,1) = sum(SigStaticClean(bin,1:31));
    PrepIndexExe(bin,1) = sum(SigStaticClean(bin,32:end));
    
    PrepExeIndex(bin,1) = sum(SigStaticClean(bin,1:31))/30;
    PrepExeIndex(bin,2) = sum(SigStaticClean(bin,32:end))/30;
    
end







