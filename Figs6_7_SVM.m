% Reproduce panels from figures 6 and 7.
% Figure 1: static bins from Grip cross-temporal decoding.
% Figure 2: cross-temporal decoding of Grip.
% Figure 3: static bins from Force cross-temporal decoding.
% Figure 4: cross-temporal decoding of Force.
% Figure 5: Decoding accuracy of grip and force (diagonals of figures 2 and
% 4).
% Figure 6: Generalization index as a function of test time.
% feb 2025, @hm

clear
close all
load DataMonkeyN
%load DataMonkeyL

NumConsecElim = 2;

Tbin = size(testvGvalidationAccuracy,1);
for bin = 1:Tbin        
       diatestvGvalidationAccuracy(bin,1) = mean(testvGvalidationAccuracy{bin,bin});
       diatestvFvalidationAccuracy(bin,1) = mean(testvFvalidationAccuracy{bin,bin});
    
end

[SigStaticGClean,testvGvalidationAccuracyClean] = CleanSigMatrix(SigStaticG,testvGvalidationAccuracy,NumConsecElim);
[StaticIndexG,SigStaticIndexG,PrepExeIndexG] = ComputeStaticIndex(SigStaticGClean,bootGvalidationAccuracy,testvGvalidationAccuracy);

[SigStaticFClean,testvFvalidationAccuracyClean] = CleanSigMatrix(SigStaticF,testvFvalidationAccuracy,NumConsecElim);
[StaticIndexF,SigStaticIndexF,PrepExeIndexF] = ComputeStaticIndex(SigStaticFClean,bootFvalidationAccuracy,testvFvalidationAccuracy);



figure
ttemcluster = Gclusters{1};
tclustersG = zeros(Tbin,1);
tclustersG(ttemcluster) = 1;

ttemcluster = Fclusters{1};
tclustersF = zeros(Tbin,1);
tclustersF(ttemcluster) = 1;
hold on
a1 = plot(xx,diatestvGvalidationAccuracy,'r-'); hold on
a2 = plot(xx,diatestvFvalidationAccuracy,'b-'); hold on
plot(xx(logical(tclustersG)),1,'r*'); hold on
plot(xx(logical(tclustersF)),1,'b*'); hold on
legend([a1 a2],{'Grip','Force'},'Location','northeast')
xlabel('Time bin #'), ylabel('Decoding accuracy')
axis tight

figure
xx = 1:1:Tbin;
hold on
plot(xx,PrepExeIndexG(:,1),'b-','linewidth',1.5);
plot(xx,PrepExeIndexG(:,2),'r-','linewidth',1.5);
axis([1, 61,0,0.61])
set(gca, 'TickDir', 'out');                     %  switch side of axis for tick marks
set(gca, 'TickLength', [.01 .01],'linewidth',1.3);            %  this in proportion of x axis
set(gca, 'XTick', [10:10:60]);        %  reset tick spacing
set(gca, 'YTick', [0:0.2:0.8]);        %  reset tick spacing
xlabel('Test time bin #'), ylabel('Generalization index')
legend({'Trained during preparation','Trained during execution'},'Location','northwest')

hold off


