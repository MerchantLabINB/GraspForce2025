function [raster] = getRasterPlots(filename)

% Plot raster plot of desired neuron.

% A list with neuron IDs for each monkey is provided in the .mat file named
% 'spikeNames.mat'.
%
% INPUTS
% filename: filename of .mat file containing spike and event times.
% two files: 'i140703-001.mat' (monkey N)
%            'l101210-001.mat' (monkey L)
% 
% spike: spike ID. ID's of example neurons shown in the paper are: 
%        Monkey N: 'spike27_1' and 'spike20_1'
%        Monkey L: 'spike27_1' and 'spike30_1'
%
% OUTPUTS
% raster: figure that contains the raster plot (top panel) and the mean
% firing rate per condition of desired neuron.
% 
% EXAMPLE
% Plot spike 20_1 of monkey N
% [raster] = getRasterPlots('i140703-001_spike27_1.mat)
% 
% feb2025, @apms. 

%  Load file
load(filename,'alignedData')

gripTypes  = {'SG','PG'};
forceTypes = {'LF','HF'};
spike = filename(13:21);

% Plot rasterplot
% Trial types are sorted from bottom to top as: 
% side grip, low force (sglf)
% side grip, high force (sghf)
% precision grip, low force (pglf)
% precision grip, high force (pghf)


conditionColor = [0.85 0 0; 0.9000 0.6000 0];
events = {'movCue','objectTouch','dispOn'};
eventsColors = [0.5 0 0; 1 0 0; 0.9, 0.5, 1];

disp_y = 0;
disp_x = 0;


% axes properties
togglefig('raster'); clf

rasterGraph = subplot(211);      % raster subplot handle
frateGraph  = subplot(212);      % firing rate subplot handle

% set raster plot parameters
rasterGraph.Position = [0.105 0.39 0.85 0.55];
rasterGraph.XLim = [-1 1];
rasterGraph.YLim = [0 50];
rasterGraph.XLimMode = 'manual';
rasterGraph.YLimMode = 'manual';
rasterGraph.XTick = [];
rasterGraph.XTickMode = 'manual';
rasterGraph.YTick = [];
rasterGraph.YTickMode = 'manual';
set(rasterGraph,'visible','off');
rasterGraph.NextPlot = 'add'; 

frateGraph.Position = [0.105 0.11 0.85 0.25];
frateGraph.NextPlot = 'add';
frateGraph.XLim = [-1 1];
frateGraph.YLim = [0 50];
frateGraph.XLimMode = 'manual';
frateGraph.YLimMode = 'manual';
frateGraph.XTick = [];
frateGraph.XTickMode = 'manual';

trialEnd = 1;

for g = 1:length(gripTypes)
    grip = gripTypes{g};

    for f = 1:length(forceTypes)
        force = forceTypes{f};

        % separate trials according to task condition
        trialTypeIndexes = strcmpi({alignedData.trials.trialType},strcat(grip,force));
        trialsToGraph = alignedData.trials(trialTypeIndexes);

        start_YTicks = 1:length(trialsToGraph);

        % obtain spike times of current task condition
        spikeTimes = alignedData.spikes.(spike);
        spikeTimes = spikeTimes(trialTypeIndexes);

        % get firing rates
        timeSamples = -1: 0.02 :trialEnd;
        fRates = firingrate(spikeTimes,timeSamples,'FilterType','exponential','TimeConstant',0.1);


        % get ticks for raster plot
        [xTicks,yTicks] = rasterplot(spikeTimes,'xlim',[-1,trialEnd],'QuickPlot','yes','displace',disp_x);

        % plot raster
        raster = togglefig('raster');
        plot(rasterGraph,xTicks,yTicks + disp_y,'.','color','k','markersize',3.5), hold on

        % plot behavioral event markers
        for m = 1:length(events)
            eventToPlot = events{m};
            plot(rasterGraph,[trialsToGraph.(eventToPlot)],(1:length(trialsToGraph))...
                +disp_y,'.','color',eventsColors(m,:),'markersize',8), hold on
        end


        % plot firing rates
        if f == 1
            lineStyle = '-';
        else
            lineStyle = '--';
        end

        plot(frateGraph,timeSamples+disp_x, mean(fRates,'omitnan'), 'color',...
            conditionColor(g,:),'linewidth',1.5,'linestyle',lineStyle), hold on

        yTicks_accum = start_YTicks + disp_y;
        disp_y = max(yTicks_accum)+8;



    end
end

% equal axes for both subplots
axis(rasterGraph,'tight');
axis(frateGraph,'tight');
rasterGraph.XLim = frateGraph.XLim;

% set xticks of firing rate graph
frateGraph.XTick = -1:0.5:1;

% set other properties for both graphs
plot(rasterGraph,[0 0],rasterGraph.YLim,'color',[0.35 0.35 0.35],'linewidth',1.5), hold on
xlabel('Time (s)'), ylabel('FR (Hz)')
plot(frateGraph,[0 0],ylim,'color',[0.35 0.35 0.35],'linewidth',1.5), hold on
legend(frateGraph,{'SGLF','SGHF','PGLF','PGHF'},'location','southwest')

rasterName = strcat(filename(1:7),'-',spike);
title(rasterGraph,rasterName);


end






