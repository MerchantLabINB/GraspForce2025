%% Plot heatmaps of auROC (of significant neurons only)

clear,clc,close all
load('ROC_data')
load('eventTimes.mat');
 
monkeys = {'monkeyN','monkeyL'};
epochEnd = 2.5;
normalize = 1;
variables = {'grip','force'};

time = -3: 0.02: 3;
xTicksTimes = -2:2;
[~,xTicks] = ismember(xTicksTimes,time);

binZero = find(time == 0);
sp1 = 0;
sp2 = 4;
sp3 = 0;

dist_neuronPerm = [];

for m = 1:length(monkeys)
    monkey = monkeys{m};
    
    if strcmpi(monkey,'monkeyN')    
        % get time bins where the mean of events occured.
        eventTimes        = eventMeanTimes.MonkeyN;
        [~,gripCueOn]     = min(abs(time - eventTimes(3)));
        [~,gripCueOff]    = min(abs(time - eventTimes(4)));
        [~,dispOn]        = min(abs(time - eventTimes(8)));
        % [~,reward]        = min(abs(time - eventTimes(9)));
        [~,objectRelease] = min(abs(time - eventTimes(10)));
    else    
        eventTimes        = eventMeanTimes.MonkeyL;
        [~,gripCueOn]     = min(abs(time - eventTimes(3)));
        [~,gripCueOff]    = min(abs(time - eventTimes(4)));
        [~,dispOn]        = min(abs(time - eventTimes(8)));
        % [~,reward]        = min(abs(time - eventTimes(9)));
        [~,objectRelease] = min(abs(time - eventTimes(10)));
    end
         
        
        if strcmpi(monkey,'monkeyN')
            meanColor = [1 0.05 0];
        else 
            meanColor = [0 0.5 1];  
        end
                  

        % Index of last time bin
        lastBin = find(round(time,2) == epochEnd);
             
        
        % ROC index     
        for v = 1:length(variables)
                variable = variables{v};
                ROC = allNeuronsROC.(monkey).(variable).ROC;              
                
                ROC = ROC(:,1:lastBin);
                
                % select neurons with significant coding only
                signifNeurons = logical(allNeuronsROC.(monkey).(variable).isSignif);
                ROC_significant = ROC(signifNeurons,:);


                % Significant neurons
                dataForTrialSorting = abs(ROC_significant-0.5); % in movIni
                [~,max_idx] = max(dataForTrialSorting,[],2);
                [~,idx] = sort(max_idx);
                sorted_ROC = dataForTrialSorting(idx,:);
                
                                              
                % Plot in the corresponding figure
                togglefig('roc heatmaps')
                
                if strcmpi(monkey,'monkeyN') && strcmpi(variable,'grip')                   
                    yLim = [0.02,0.25];
                    clims = [0 0.4];
                                                    
                elseif strcmpi(monkey,'monkeyN') && strcmpi(variable,'force')
                    yLim = [0.02,0.25];
                    clims = [0 0.4];

                    
                elseif strcmpi(monkey,'monkeyL') && strcmpi(variable,'grip')
                    yLim = [0.02,0.25];
                    clims = [0 0.4];
               
                elseif strcmpi(monkey,'monkeyL') && strcmpi(variable,'force')
                    yLim = [0.02,0.25];
                    clims = [0 0.4];

                end
                
                
                % Plot matrix of ROC of all neurons
                subplot(2,4,sp1+v)               
                maxBin = 76;
                imagesc([sorted_ROC(:,1:lastBin)],clims), hold on
                line([binZero, binZero],ylim,'color','k','linewidth',1)
                line([gripCueOn, gripCueOn],ylim,'color',[0.6 0.6 0.6],'linewidth',1), hold on
                line([gripCueOff, gripCueOff],ylim,'color',[0.6 0.6 0.6],'linewidth',1), hold on
                line([dispOn, dispOn],ylim,'color',[0.6 0.6 0.6],'linewidth',1), hold on
                line([objectRelease, objectRelease],ylim,'color',[0.6 0.6 0.6],'linewidth',1)
                xticks(xTicks);
                xticklabels({xTicksTimes})
                title(strcat("Monkey ",monkey(end)," - ",variable))
                ylabel('Cell #'), box off, xlim([51 251])
                % colorbar
                
                % Plot mean of ROC                
                dataForPlotting = ROC_significant(:,1:lastBin); 
                mean_ROC = mean(abs(dataForPlotting-0.5));
                sem_ROC = std(abs(dataForPlotting-0.5))./sqrt(size(dataForPlotting,1));
                
                subplot(2,4,sp2+v)
                plot(time(1:lastBin),mean_ROC,'color',meanColor,'linewidth',1.3), hold on
                myerrorbar(time(1:lastBin),mean_ROC,sem_ROC,...
                    'color',meanColor,'interval',1,'alpha',0.25), hold on
                
                xlim([-2 2]), ylim(yLim)
                line([0 0],ylim,'color','k','linewidth',1); hold on   
                line([eventTimes(8), eventTimes(8)],ylim,'color',[0.6 0.6 0.6],'linewidth',1), hold on
                line([eventTimes(10), eventTimes(10)],ylim,'color',[0.6 0.6 0.6],'linewidth',1)
                ylabel('abs(auROC)'), box off 

                                                                                  
                clear pointsAboveCi
                clear significantBins
                clear shuff_codingProportion
                clear codingProportion
                signifNeurons = [];
                activationPeriod = [];
                latency = [];
        end


    
    sp1 = 2;
    sp2 = 6;
    sp3 = 1;
    
end
