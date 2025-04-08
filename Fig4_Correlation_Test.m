% Perform non parametric correlation test in neurons with significant auROC 
% (grip-force and grip-only neurons).
% feb 2025, @apms
%

clear,clc,close all
load('correlationPermutations.mat')
load('ROC_data.mat')
load('eventTimes')

monkeys = {'monkeyN','monkeyL'};
timeAxis = -3:0.02:3;

for m = 1:length(monkeys)
    monkey = monkeys{m};

    if m == 1
        ot_timebin = find(round(timeAxis,2) == 0.4); % object touch bin
        ot_time = eventMeanTimes.MonkeyN(7);
    else
        ot_timebin = find(round(timeAxis,2) == 0.16); % object touch bin
        ot_time = eventMeanTimes.MonkeyL(7);
    end

    grip_units = logical(allNeuronsROC.(monkey).grip.isSignif);
    force_units = logical(allNeuronsROC.(monkey).force.isSignif);

    % indexes of only-grip and grip-force selective units
    onlygrip_units = grip_units & ~force_units;
    gripforce_units = grip_units & force_units;

    % auROC values
    grip_ROC = allNeuronsROC.(monkey).grip.ROC;
    force_ROC = allNeuronsROC.(monkey).force.ROC;

    % absolute auROC values
    abs_grip_ROC = abs(allNeuronsROC.(monkey).grip.ROC - 0.5); 
    abs_force_ROC = abs(allNeuronsROC.(monkey).force.ROC - 0.5);

   
    % Plot abs(auROC) as a function of time (significant neurons)

    % ---------------- Only-grip neurons ----------------------------------
    togglefig(monkey)
    subplot(2,3,[1 2])

    % plot grip type auROC
    GripROC_mean = mean(abs_grip_ROC(onlygrip_units,:),'omitnan');
    GripROC_stdError = std(abs_grip_ROC(onlygrip_units,:),[],'omitnan')/...
        sqrt(size(abs_grip_ROC,1));

    plot(timeAxis, GripROC_mean,'-r', 'LineWidth',1.5), hold on
    myerrorbar(timeAxis, GripROC_mean, GripROC_stdError,'color','r',...
        'interval',1), hold on
    
    % plot force level auROC
    ForceROC_mean = mean(abs_force_ROC(onlygrip_units,:),'omitnan');
    ForceROC_stdError = std(abs_force_ROC(onlygrip_units,:),[],'omitnan')...
        /sqrt(size(abs_force_ROC,1));
    
    plot(timeAxis, ForceROC_mean,'-b','LineWidth',1.5), hold on
    myerrorbar(timeAxis, ForceROC_mean, ForceROC_stdError,'color','b',...
        'interval',1), hold on    

    % mark significantly dependent bins
    onlygrip_dependent_bins = logical(correlationPermutations.(monkey).onlygrip_units.correlated_bins);
    if sum(onlygrip_dependent_bins) > 0
        plot(timeAxis(onlygrip_dependent_bins),0.25,'*k'), hold on
    end

    title(strcat("only grip neurons ",num2str(sum(onlygrip_units)),'/',...
        num2str(length(onlygrip_units)))), xlim([-2 2]);

    plot([ot_time ot_time],ylim,'--k'), hold on
    xlabel('Time (s)'), ylabel('abs(auROC)'), axis tight
    
   
    % Plot correlation coefficients and linear regression
    subplot(233)

    x = abs_force_ROC(onlygrip_units,ot_timebin);
    y = abs_grip_ROC(onlygrip_units,ot_timebin);

    % fit simple linear regression
    p = polyfit(x,y,1);
    f = polyval(p,x);

    onlygrip_pvalues = correlationPermutations.(monkey).onlygrip_units.pValues;
    onlygrip_corrcoefs = round(correlationPermutations.(monkey).onlygrip_units.rho,2);

    plot(x, y, '.', x,f,'-'), hold on
    text(0.15,0.4,strcat("p = ",num2str(onlygrip_pvalues(ot_timebin)))), hold on
    text(0.15,0.2,strcat("R = ",num2str(onlygrip_corrcoefs(ot_timebin)))), hold on
    xticks(0:0.1:0.5), yticks(0:0.1:0.5)
    title(strcat("OT time ",'(',num2str(timeAxis(ot_timebin)*1000)," ms",')'))
    xlabel('abs auROC (force)'), ylabel('abs auROC (grip)');


    

    % ----------------------- Grip-force neurons --------------------------
    subplot(2,3,[4 5])
    
    % plot grip type auROC
    GripROC_mean = mean(abs_grip_ROC(gripforce_units,:),'omitnan');
    GripROC_stdError = std(abs_grip_ROC(gripforce_units,:),[],'omitnan')/...
        sqrt(size(abs_grip_ROC,1));

    plot(timeAxis, GripROC_mean,'-r', 'LineWidth',1.5), hold on
    myerrorbar(timeAxis, GripROC_mean, GripROC_stdError,'color','r',...
        'interval',1), hold on
    
    % plot force level auROC
    ForceROC_mean = mean(abs_force_ROC(gripforce_units,:),'omitnan');
    ForceROC_stdError = std(abs_force_ROC(gripforce_units,:),[],'omitnan')...
        /sqrt(size(abs_force_ROC,1));
    
    plot(timeAxis, ForceROC_mean,'-b','LineWidth',1.5), hold on
    myerrorbar(timeAxis, ForceROC_mean, ForceROC_stdError,'color','b',...
        'interval',1), hold on
      

    % mark significantly dependent bins
    gripforce_dependent_bins = logical(correlationPermutations.(monkey).gripforce_units.correlated_bins);
    if sum(gripforce_dependent_bins) > 0
        plot(timeAxis(gripforce_dependent_bins),0.25,'*k'), hold on
    end
    
    title(strcat("grip/force neurons ",num2str(sum(gripforce_units)),'/',...
        num2str(length(gripforce_units)))), xlim([-2 2]);
   
    plot([ot_time ot_time],ylim,'--k'), hold on
    xlabel('Time (s)'), ylabel('abs(auROC)'), axis tight

   
    % Plot correlation coefficients and linear regression
    subplot(236)
    x = abs_force_ROC(gripforce_units,ot_timebin);
    y = abs_grip_ROC(gripforce_units,ot_timebin);

    % fit simple linear regression
    p = polyfit(x,y,1);
    f = polyval(p,x);

    gripforce_pvalues = correlationPermutations.(monkey).gripforce_units.pValues;
    gripforce_corrcoefs = round(correlationPermutations.(monkey).gripforce_units.rho,2);

    plot(x, y, '.', x,f,'-'), hold on
    text(0.15,0.4,strcat("p = ",num2str(gripforce_pvalues(ot_timebin)))), hold on
    text(0.15,0.2,strcat("R = ",num2str(gripforce_corrcoefs(ot_timebin)))), hold on
    xticks(0:0.1:0.5), yticks(0:0.1:0.5)
    title(strcat("OT time ",'(',num2str(timeAxis(ot_timebin)*1000)," ms",')'))
    xlabel('abs auROC (force)'), ylabel('abs auROC (grip)');



end



    

        