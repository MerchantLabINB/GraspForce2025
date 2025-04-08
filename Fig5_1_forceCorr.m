
% Figure 5-1.
% Estimate the correlation between condition-independent and force-related 
% dPCs and the pulling force signal (trial by trial).
% @apms, 2025

monkeys = {'monkeyN','monkeyL'};
load('force_signals')
FR_timeAxis = -3:0.02:3; % time axis of firing rates

conditions = {'SGLF','SGHF','PGLF','PGHF'};
conditionColors = [0 0.2 0.9; 0 0.2 0.9; 0 0.8 0.5; 0 0.8 0.5];
forceColors = [0 0 1; 1 0 0];

allTrialsForce = [];
allTrialsDisp = [];

xcf_grip = nan(45,4);
xcf_force = nan(45,4);
xcf_time = nan(45,4);


dpca_margs = {'force','time'}; 

cont = 0;

for m = 1:length(monkeys)
    monkey = monkeys{m};

    if m == 1 
        load('dPCA_weighted_averages_monkeyN.mat')
    else
        load('dPCA_weighted_averages_monkeyL.mat')
    end

    data = eval(monkey);

    % select hits only
    hits = ~cellfun(@isempty,{data.trials.objectTouch});

    trials = data.trials(hits);
    trials_force = data.force_signals.pullingForce(hits);
    trials_disp = data.force_signals.objectDisp(hits);
    
    Fs = 0.02; % downsampled
    meanTimeAxis = -2:Fs:200*Fs;
    meanTimeAxis = meanTimeAxis(1:end-1);

    % Mean times of other important events
    objectTouch = mean([trials.objectTouch]-[trials.movIni]);
    dispOn = mean([trials.dispOn] - [trials.movIni]);
    objectRelease = mean([trials.objectRelease]-[trials.movIni]);

    % separate trials per condition
    for c = 1:length(conditions)
        trialType = conditions{c};
        trialTypeIdx = ismember({trials.trialType},trialType) ;

        conditionTrials = trials(trialTypeIdx);

        conditionForce = trials_force(trialTypeIdx);
        conditionDisp = trials_disp(trialTypeIdx); 

        fillUntil = max(cellfun(@length,conditionForce));

        % Attrition 
        conditionForce = cell2mat(cellfun(@(x) fillWithNans(x, fillUntil), conditionForce, 'UniformOutput',false));
        conditionDisp = cell2mat(cellfun(@(x) fillWithNans(x, fillUntil), conditionDisp, 'UniformOutput',false));

        % Time axis
        timeAxis = 0:Fs:length(conditionForce)*Fs;
        timeAxis = timeAxis(1:end-1);

        % Align trials to movIni
        movIni = [conditionTrials.movIni];
        trialStarts = movIni - 2;

        conditionForceAligned = nan(length(movIni),201);
        conditionDispAligned = nan(length(movIni),201);

        for t = 1:length(movIni)
            trialBins = timeAxis >= trialStarts(t);

            % align force signal
            trialForceAligned = conditionForce(t,trialBins);
            
            if length(trialForceAligned) < 201
                trialForceAligned(end+1:201) = nan;
            else
                trialForceAligned = trialForceAligned(1:201);
            end


            % align displacement signal
            trialDispAligned = conditionDisp(t,trialBins);
            
            if length(trialDispAligned) < 201
                trialDispAligned(end+1:201) = nan;
            else
                trialDispAligned = trialDispAligned(1:201);
            end
            
            conditionForceAligned(t,:) = trialForceAligned; 
            conditionDispAligned(t,:) = trialDispAligned; 
        end

        alignedTimeAxis = -2:0.02:2;

        for mg = 1:length(dpca_margs)
            marg = dpca_margs{mg};

            weightedMeans = weightedAverages.(trialType).(marg);

            for t = 1:size(weightedMeans,1)
                trialWeightedAverage = weightedMeans(t,:);
                
                % make weighted average the same length as force signal
                trialWeightedAverage = trialWeightedAverage(51:251);

                trialForceSignal = conditionForceAligned(t,:);

                % % plot for inspection
                % figure
                % yyaxis left
                % plot(alignedTimeAxis,trialWeightedAverage), hold on
                % yyaxis right
                % plot(alignedTimeAxis,trialForceSignal), hold on

                % estimate trial cross-correlation
                binstokeep = ~isnan(trialForceSignal);

                [xcf,lags,bounds] = crosscorr(trialWeightedAverage(binstokeep),trialForceSignal(binstokeep),NumSTD=3,NumLags=10);
                signifCoeffs = xcf > bounds(1) | xcf < bounds(2); 

                crossCorrelation.xcf(t,:) = abs(xcf);
                crossCorrelation.lags(t,:) = lags;
                crossCorrelation.isSignif(t,:) = signifCoeffs;


            end


            signifCorrProportion = sum(crossCorrelation.isSignif,1)/size(crossCorrelation.isSignif,1);

            [max_xcf,max_idx] = max(crossCorrelation.xcf,[],2);
            max_xcf_lag = lags(max_idx);
            max_signif_trials = double(diag(crossCorrelation.isSignif(:,max_idx)));
            max_signif_trials(end+1:45) = nan;

            
            if c == 2 || c == 4
                lnstyle = '--';    
            else
                lnstyle = '-';
            end


            if mg == 1
                xcf_force(1:length(max_xcf),c) = abs(max_xcf);
                 prop_force(:,c) = max_signif_trials;
            elseif mg == 2
                 xcf_time(1:length(max_xcf),c) = abs(max_xcf);
                 prop_time(:,c) = max_signif_trials;
            end
            
            
            % normalize force signal to zscore
            zmean = mean(conditionForceAligned(:,1:101),2);
            zstd = std(conditionForceAligned(:,1:101),0,2);
            norm_conditionForceAligned = (conditionForceAligned - zmean)./zstd;

            timeAxis2 = -2:0.02:2;
            togglefig(strcat(monkey,marg))
            subplot(221)          
            plot(timeAxis2,mean(norm_conditionForceAligned,'omitnan'),'color',conditionColors(c,:),'LineStyle',lnstyle,'LineWidth',1.5), hold on
            ylabel('Force (zscore)'), box off, title('Pulling-force signal')

            subplot(222)
            plot(timeAxis2,mean(weightedMeans(:,51:251)),'color',conditionColors(c,:),'LineStyle',lnstyle,'LineWidth',1.5), hold on
            ylabel('Norm. firing rate (Hz)'), box off, title(strcat(marg," dPC"))

            subplot(223)
            plot(-10:10,signifCorrProportion,'color',conditionColors(c,:),'LineStyle',lnstyle,'LineWidth',1.5), hold on
            xlabel('Lags'), ylabel('Proportion of correlated trials'), box off
            title(marg)

            subplot(224)
            plot(-10:10, mean(crossCorrelation.xcf),'color',conditionColors(c,:),'LineStyle',lnstyle,'LineWidth',1.5), hold on
            xlabel('Lags'), ylabel('Cross-correlation coefficient'), box off

            legend({'Side grip, LF','Side grip, HF', 'Prec. grip LF','Prec. grip, HF'},'location','northwest')

          
        end
        
        cont = cont + 5;
        conditionForceAligned = [];
        conditionDispAligned = [];
        
        
    end

    % perform anova
    
    xcf_force_low = [rmmissing(xcf_force(:,1)); rmmissing(xcf_force(:,3))];
    xcf_force_high = [rmmissing(xcf_force(:,2)); rmmissing(xcf_force(:,4))];

    xcf_time_low = [rmmissing(xcf_time(:,1)); rmmissing(xcf_time(:,3))];
    xcf_time_high = [rmmissing(xcf_time(:,2)); rmmissing(xcf_time(:,4))];


    anovaVector = [xcf_force_low; xcf_force_high; xcf_time_low; xcf_time_high];

    forceLevel = [ones(length(xcf_force_low),1); ones(length(xcf_force_high),1)*2; ones(length(xcf_time_low),1); ones(length(xcf_time_high),1)*2];

    marginalization = [ones(length(xcf_force_low),1)*2;...
        ones(length(xcf_force_high),1)*2; ones(length(xcf_time_low),1)*3; ones(length(xcf_time_high),1)*3];

    
    figure
    [~,~,stats] = anovan(anovaVector,{forceLevel,marginalization},"model","interaction","varnames",["forceLevel","marginalization"]);
    results = multcompare(stats,"Dimension",2);


    

    
    prop_force_LF = [rmmissing(prop_force(:,1)); rmmissing(prop_force(:,3))];
    prop_force_LF = sum(prop_force_LF)/length(prop_force_LF);

    prop_force_HF = [rmmissing(prop_force(:,2)); rmmissing(prop_force(:,4))];
    prop_force_HF = sum(prop_force_HF)/length(prop_force_HF);

    
    prop_time_LF = [rmmissing(prop_time(:,1)); rmmissing(prop_time(:,3))];
    prop_time_LF = sum(prop_time_LF)/length(prop_time_LF);

    prop_time_HF = [rmmissing(prop_time(:,2)); rmmissing(prop_time(:,4))];
    prop_time_HF = sum(prop_time_HF)/length(prop_time_HF);

    togglefig(strcat(monkey,'-boxplots'))
    subplot(122)
    y = [prop_force_LF prop_force_HF; prop_time_LF prop_time_HF];
    bar(y), ylabel('proportion of trials'), box off


    xcf_margs = {{'xcf_force_low','xcf_force_high'},{'xcf_time_low','xcf_time_high'}};

    for n = 1:length(xcf_margs)
        xcf_low = eval(xcf_margs{1,n}{1,1});
        xcf_high = eval(xcf_margs{1,n}{1,2});

        subplot(121)
        group = repmat(n+cont,length(xcf_low),1);
        b = boxchart(group, xcf_low); hold on
        b.JitterOutliers = 'on';
        b.MarkerStyle = '.';
        b.BoxFaceColor = forceColors(1,:);

        group = repmat(n+cont+1,length(xcf_high),1);
        b = boxchart(group, xcf_high); hold on
        b.JitterOutliers = 'on';
        b.MarkerStyle = '.';
        b.BoxFaceColor = forceColors(2,:);

        if n == 1
            marg = 'force';
        elseif n == 2
            marg = 'time';
        end

        text(n+cont+1-0.5,1,marg)
        cont = cont+4; 
    end
  
    xticks([])
    ylabel('absolute xcf')
    ylim([0 1])
    title(monkey)
    legend({'low force', 'high force'},'Location','southeast')

end



function c = fillWithNans(c,limit)
c(end:limit) = nan;
end

