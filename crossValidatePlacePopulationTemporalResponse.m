function crossValidatedStats = crossValidatePlacePopulationTemporalResponse(firing_rates,num_shuffs,type,twin1,twin2)
%turned into function on  5/21/20
%firing rates are cell array of smoothed firing rates by trial/fixation
%type is place vs nonplace


%store shuffled firing rates so don't have to do this more than once
num_units = length(firing_rates);


%---Get Temporal Stability Measurements For Each Place Cell & Shuffle Firing Rates---%
observed_distance = NaN(1,num_units); %distance between peaks
observed_COM = NaN(1,num_units);%distance between coms
observed_corr = NaN(1,num_units); %correlation between firing rate curves
shuff_distance = NaN(num_shuffs,num_units); %shuffled peak distances
shuff_COM = NaN(num_shuffs,num_units); %shuffled COM distances
shuff_corr = NaN(num_shuffs,num_units); %shuffled firing rate corrs

all_shuff_firing_rates = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates{i,shuff} = NaN(num_units,twin1+twin2);
    end
end

for unit = 1:num_units
    
    Dmean1 = mean(firing_rates{unit}(1:2:end,:));
    Dmean2 = mean(firing_rates{unit}(2:2:end,:));
    
    %center of mass firing rate
    com1 = sum(Dmean1.*[1:twin1+twin2])/sum(Dmean1);
    com2 = sum(Dmean2.*[1:twin1+twin2])/sum(Dmean2);
    
    %get peak times, not using other more complicated methods
    [~,mx1] = max(Dmean1,[],2);
    [~,mx2] = max(Dmean2,[],2);
    
    observed_COM(unit) = abs(com1-com2);
    observed_distance(unit) = abs(mx1-mx2); 
    observed_corr(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        %shuffle firing times, row shuffling creates flat firing rate with
        %center of mass around center of time window
        Dmean1 = circshift_row(Dmean1);
        Dmean2 = circshift_row(Dmean2);
        com1 = sum(Dmean1.*[1:twin1+twin2])/sum(Dmean1);
        com2 = sum(Dmean2.*[1:twin1+twin2])/sum(Dmean2);
        
        %shuffle firing rates
        shuff_firing = circshift_row(firing_rates{unit});
        
        Dmean11 = mean(shuff_firing(1:2:end,:));
        Dmean22 = mean(shuff_firing(2:2:end,:));

        [~,mx1] = max(Dmean11,[],2);
        [~,mx2] = max(Dmean22,[],2);

        shuff_COM(shuff,unit) = abs(com1-com2);
        shuff_distance(shuff,unit) = abs(mx1-mx2); 
        shuff_corr(shuff,unit) = corr(Dmean11',Dmean22','row','pairwise','type','Spearman');
        
        all_shuff_firing_rates{1,shuff}(unit,:) = Dmean11;
        all_shuff_firing_rates{2,shuff}(unit,:) = Dmean22;
    end
end

%get proportion of significant units
prctile_distance = NaN(1,num_units);
prctile_com = NaN(1,num_units);
prctile_corr = NaN(1,num_units);
for unit = 1:num_units
    prctile_corr(unit) = 100*sum(observed_corr(unit) > shuff_corr(:,unit))/num_shuffs;
    prctile_distance(unit) = 100*sum(observed_distance(unit) < shuff_distance(:,unit))/num_shuffs; 
    prctile_com(unit) = 100*sum(observed_COM(unit) < shuff_COM(:,unit))/num_shuffs; 
end


%---Get Population-Level Correlation---%
firing_rates1 = NaN(num_units,twin1+twin2);
firing_rates2 = NaN(num_units,twin1+twin2);
centerMass1 = NaN(num_units,1);
centerMass2 = NaN(num_units,1);
for unit = 1:num_units
    
    Dmean1 = mean(firing_rates{unit}(1:2:end,:));
    Dmean2 = mean(firing_rates{unit}(2:2:end,:));
    
    centerMass1(unit) = sum(Dmean1.*[1:twin1+twin2])/sum(Dmean1);
    centerMass2(unit) = sum(Dmean2.*[1:twin1+twin2])/sum(Dmean2);
    
    Dmean1 = Dmean1-mean(Dmean1(1:twin1));
    Dmean1 = Dmean1/max(Dmean1);
    
    Dmean2 = Dmean2-mean(Dmean2(1:twin1));
    Dmean2 = Dmean2/max(Dmean2);
    
    firing_rates1(unit,:) = Dmean1;
    firing_rates2(unit,:) = Dmean2;
end

[~,mx1] = max(firing_rates1,[],2);
[~,mx2] = max(firing_rates2,[],2);

observed_peak_corrs = corr(mx1,mx2,'row','pairwise','type','Spearman');
observed_COM_corrs = corr(centerMass1,centerMass2,'row','pairwise','type','Spearman');
observed_corr_timecellplot = corr(firing_rates1(:),firing_rates2(:),'row','pairwise','type','Spearman');
observed_order_distance = nanmedian(abs(mx1-mx2));
observed_com_distance = nanmedian(abs(centerMass1-centerMass2));

%---Get Distribution of Shuffled Population-Level Correlations---%

shuff_peak_corrs =  NaN(1,num_shuffs);
shuff_corr_timecellplot =  NaN(1,num_shuffs);
shuff_order_distance = NaN(1,num_shuffs);
shuff_order_com = NaN(1,num_shuffs);
for shuff = 1:num_shuffs
    
    shuff_firining_rates1 = all_shuff_firing_rates{1,shuff};
    shuff_firining_rates2 = all_shuff_firing_rates{2,shuff};
    
    shuff_coms = NaN(2,num_units);
    
    for unit = 1:num_units
        fr = shuff_firining_rates1(unit,:);
        fr = fr-mean(fr(1:twin1));
        fr = fr/max(fr);
        shuff_firining_rates1(unit,:) = fr;
        shuff_coms(1,unit) = sum(fr.*[1:twin1+twin2])/sum(fr);
        
        fr = shuff_firining_rates2(unit,:);
        fr = fr-mean(fr(1:twin1));
        fr = fr/max(fr);
        shuff_firining_rates2(unit,:) = fr;
        shuff_coms(2,unit) = sum(fr.*[1:twin1+twin2])/sum(fr);
    end
    
    [~,shuff_mx1] = max(shuff_firining_rates1,[],2);    
    [~,shuff_mx2] = max(shuff_firining_rates2,[],2);
    
    shuff_peak_corrs(shuff) = corr(shuff_mx1,shuff_mx2,'row','pairwise','type','Spearman');
    shuff_corr_timecellplot(shuff) = corr(shuff_firining_rates1(:),shuff_firining_rates2(:),'row','pairwise','type','Spearman');
    shuff_order_distance(shuff) = median(abs(shuff_mx1-shuff_mx2));
    shuff_order_com(shuff) = median(abs(shuff_coms(1,:)-shuff_coms(2,:)));
end

%% 
%---Get Sum Stats---%
numSpikes = NaN(1,num_units);
numFixations = NaN(1,num_units);
numFixationsWSpikes = NaN(1,num_units);
for unit = 1:num_units
    numSpikes(unit) = sum(firing_rates{unit}(:));
    numFixations(unit) = size(firing_rates{unit},1);
    numFixationsWSpikes(unit) = sum(sum(firing_rates{unit}') > 0);
end

[rhoSpikesPeak,pSpikesPeak] = corr(numSpikes',prctile_distance','row','pairwise','type','Spearman');
[rhoFixPeak,pFixPeak] = corr(numFixations',prctile_distance','row','pairwise','type','Spearman');
[rhoFixSpikesPeak,pFixSpikesPeak] = corr(numFixationsWSpikes',prctile_distance','row','pairwise','type','Spearman');


[rhoSpikesCorr,pSpikesCorr] = corr(numSpikes',prctile_corr','row','pairwise','type','Spearman');
[rhoFixCorr,pFixCorr] = corr(numFixations',prctile_corr','row','pairwise','type','Spearman');
[rhoFixSpikesCorr,pFixSpikesCorr] = corr(numFixationsWSpikes',prctile_corr','row','pairwise','type','Spearman');
%%
%---Plot Results---%
figure
subplot(2,3,1)
plot(mx1-twin1,mx2-twin1,'.k')
xlabel('Peak Odd Trials (ms)')
ylabel('Peak Even Trials (ms)')
box off
title(['Population-level Peak Correlation (\rho) = ' num2str(observed_peak_corrs,3)])
axis equal

subplot(2,3,4)
plot(centerMass1-twin1,centerMass2-twin1,'.k')
xlabel('COM Odd Trials (ms)')
ylabel('COM Even Trials (ms)')
box off
title(['Population-level COM Correlation (\rho) = ' num2str(observed_COM_corrs,3)])
axis equal

subplot(2,3,2)
plot(log10(numSpikes),prctile_corr,'k.')
hold on
plot(log10(numSpikes),prctile_distance,'b.')
hold off
xlabel('Log_{10}(# Spikes)')
ylabel('Sig. Percentile')
title(['\rho_{spikes2corr} = ' num2str(rhoSpikesCorr,3) ', p = ' num2str(pSpikesCorr,3) ...
    ', \rho_{spikes2peak} = ' num2str(rhoSpikesPeak,3) ', p = ' num2str(pSpikesPeak,3)])
box off

subplot(2,3,3)
plot(numFixations,prctile_corr,'k.')
hold on
plot(numFixations,prctile_distance,'b.')
hold off
xlabel('# Fixations')
ylabel('Sig. Percentile')
title(['\rho_{fix2corr} = ' num2str(rhoFixCorr,3) ', p = ' num2str(pFixCorr,3) ...
    ', \rho_{fix2peak} = ' num2str(rhoFixPeak,3) ', p = ' num2str(pFixPeak,3)])
box off

subplot(2,3,6)
plot(numFixationsWSpikes,prctile_corr,'k.')
hold on
plot(numFixationsWSpikes,prctile_distance,'b.')
hold off
xlabel('# Fix w/ Spikes')
ylabel('Sig. Percentile')
title(['\rho_{fixS2corr} = ' num2str(rhoFixSpikesCorr,3) ', p = ' num2str(pFixSpikesCorr,3) ...
    ', \rho_{fixs2peak} = ' num2str(rhoFixSpikesPeak,3) ', p = ' num2str(pFixSpikesPeak,3)])
box off
%%

[~,i1] = sort(centerMass1);
[~,i2] = sort(centerMass2);


figure
negVal = -nanstd(firing_rates1(:));
subplot(2,3,1)
imagesc(-twin1:twin2-1,1:size(firing_rates1,1),firing_rates1(i1,:));
caxis([negVal 1])
title('Odd Trials Sorted by Odd')
xlabel('Time from Fixation (ms)')
ylabel('View Cell#')

negVal = -nanstd(firing_rates1(:));
subplot(2,3,2)
imagesc(-twin1:twin2-1,1:size(firing_rates1,1),firing_rates1(i2,:));
caxis([negVal 1])
title('Odd Trials Sorted by Even')
xlabel('Time from Fixation (ms)')
ylabel('View Cell#')

negVal = -nanstd(firing_rates2(:));
subplot(2,3,4)
imagesc(-twin1:twin2-1,1:size(firing_rates2,1),firing_rates2(i2,:));
caxis([negVal 1])
title('Even Trials Sorted by Even')
xlabel('Time from Fixation (ms)')
ylabel('View Cell#')

negVal = -nanstd(firing_rates2(:));
subplot(2,3,5)
imagesc(-twin1:twin2-1,1:size(firing_rates2,1),firing_rates2(i1,:));
caxis([negVal 1])
title('Even Trials Sorted by Odd')
xlabel('Time from Fixation (ms)')
ylabel('View Cell#')

subplot(2,3,[3 6])
plot([-twin1 twin2],[0 0],'k')
hold on
plot([0 0],[-0.2 1],'k')
plot(-twin1:twin2-1,nanmean(firing_rates1))
plot(-twin1:twin2-1,nanmean(firing_rates2))
hold off
ylim([-0.2 1])
xlabel('Time from Fixation (ms)')
ylabel('Normalized Firing Rate')
legend('Odd Trials','Even Trials')
title('Average Firing Rate')
box off
axis square


sgtitle(type)

%%
%---Store All Data in 1 Structure Array---%
crossValidatedStats = [];
crossValidatedStats.placeType = type;
crossValidatedStats.num_units = num_units;

crossValidatedStats.individualStats.observed_distance = observed_distance;
crossValidatedStats.individualStats.observed_COM = observed_COM;
crossValidatedStats.individualStats.obseser = observed_corr;
crossValidatedStats.individualStats.shuff_distance = shuff_distance;
crossValidatedStats.individualStats.shuff_COM = shuff_COM;
crossValidatedStats.individualStats.shuff_corr = shuff_corr;
crossValidatedStats.individualStats.prctile_distance = prctile_distance;
crossValidatedStats.individualStats.prctile_com = prctile_com;
crossValidatedStats.individualStats.prctile_corr = prctile_corr;

crossValidatedStats.populationStats.observed_peak_corrs = observed_peak_corrs;
crossValidatedStats.populationStats.observed_COM_corrs = observed_COM_corrs;
crossValidatedStats.populationStats.observed_corr_timecellplot = observed_corr_timecellplot;
crossValidatedStats.populationStats.observed_order_distance = observed_order_distance;
crossValidatedStats.populationStats.observed_com_distance = observed_com_distance;
crossValidatedStats.populationStats.shuff_peak_corrs = shuff_peak_corrs;
crossValidatedStats.populationStats.shuff_corr_timecellplot = shuff_corr_timecellplot;
crossValidatedStats.populationStats.shuff_order_distance = shuff_order_distance;


%%
%% Main Print Results
disp('------------------------------------------------------')
disp(['Cross Validation Stats for ' type ' cells (n = ' num2str(num_units) ')']);
disp('Inidividual Unit Stats:')
disp(['Firing Rate Correlation: ' num2str(sum(prctile_corr > 95)) '(' num2str(100*sum(prctile_corr > 95)/num_units,3) '%)'])
disp(['Peak Time Distance: ' num2str(sum(prctile_distance > 95)) '(' num2str(100*sum(prctile_distance > 95)/num_units,3) '%)'])
disp(['COM Time Distance: ' num2str(sum(prctile_com > 95)) '(' num2str(100*sum(prctile_com > 95)/num_units,3) '%)'])
disp(' ')

disp('Population-Level Stats:')
disp(['Population Level Correlation: ' num2str(observed_corr_timecellplot) ' (' num2str(100*sum(observed_corr_timecellplot > shuff_corr_timecellplot)/num_shuffs,3) '%)'])
disp(['Population Level Distance: ' num2str(observed_order_distance) ' ms (' num2str(100*sum(observed_order_distance < shuff_order_distance)/num_shuffs,3) '%)'])
disp(['Population Level COM: ' num2str(observed_com_distance) ' ms (' num2str(100*sum(observed_com_distance < shuff_order_com)/num_shuffs,3) '%)'])

disp(' ')

end