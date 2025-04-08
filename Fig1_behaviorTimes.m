% Get reaction and movement times and plot them using violin plots

load('behavTimes.mat')

reactionTimes{1} = B.monkeyN.reaction_times;
reactionTimes{2} = B.monkeyL.reaction_times;
movementTimes{1} = B.monkeyN.movement_times;
movementTimes{2} = B.monkeyL.movement_times;

% Plot violins using violin.m function:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de

figure,clf
subplot(121)
[h,L,MX,MED] = violin(reactionTimes,'facecolor',[1 0 0; 0 0.4 1],...
    'edgecolor','none','xlabel',{"Monkey N", "Monkey L"});
ylabel('Reaction time (s)'), grid on

subplot(122)
[h,L,MX,MED] = violin(movementTimes,'facecolor',[1 0 0; 0 0.4 1],...
    'edgecolor','none','xlabel',{"Monkey N", "Monkey L"});
ylabel('Movement time (s)'), grid on