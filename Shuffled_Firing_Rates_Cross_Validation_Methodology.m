clar
Fs = 1000;
smth = 30;

num_shuffs = 1000;

num_units = 100;
twin1 = -200;
twin2 = 400;2
maxt = twin2-twin1;

corr50 = NaN(num_shuffs,num_units);
corr_index = NaN(1,num_shuffs);
corr_peaks = NaN(1,num_shuffs);
for shuff = 1:num_shuffs
    all_firing_rates = cell(1,num_units);
    all_mean_firing_rates = NaN(num_units,maxt);
    
    all_mean_firing_rates1 = NaN(num_units,maxt);
    all_mean_firing_rates2 = NaN(num_units,maxt);

    for i = 1:num_units
        
        num_trials = randi([100 2000]);
        z = randn(num_trials,maxt);
        
        [Dmean, Dtrls] = nandens (z, smth, 'gauss', Fs);
        all_firing_rates{i} = Dtrls;
        all_mean_firing_rates(i,:) = Dmean;
        
%         Dmean1 = mean(Dtrls(1:floor(num_trials/2),:));
%         Dmean2 = mean(Dtrls(ceil(num_trials/2):end,:));
        Dmean1 = mean(Dtrls(1:2:end,:));
        Dmean2 = mean(Dtrls(2:2:end,:));

        
        all_mean_firing_rates1(i,:) = Dmean1;
        all_mean_firing_rates2(i,:) = Dmean2;
        
        corr50(shuff,i) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    end    
    
    [~,mx1] = max(all_mean_firing_rates1,[],2);
    [~,ii1] = sort(mx1);
    
    [~,mx2] = max(all_mean_firing_rates2,[],2);
    [~,ii2] = sort(mx2);
    
    corr_index(shuff) = corr(ii1,ii2,'row','pairwise','type','Spearman');
    corr_peaks(shuff) = corr(mx1,mx2,'row','pairwise','type','Spearman');
end
%%
% figure
% subplot(2,2,1)
% imagesc(twin1:twin2-1,1:num_units,all_mean_firing_rates);
% title('Smoothed-Raw')
% axis square
% colorbar
% caxis([pct25 pct975])
% 
% subplot(2,2,2)
% imagesc(twin1:twin2-1,1:num_units,all_mean_firing_rates(ii,:));
% title('Sorted-Smoothed Raw')
% colormap('jet')
% axis square
% colorbar
% caxis([pct25 pct975])
% 

% subplot(2,2,3)
% hist(corr50(:),50)
% box off
% xlabel('Corr_{1/2}')
% ylabel('# Units')
% title(['Medan Corr_{1/2} = ' num2str(mean(corr50),3)])

% subplot(2,2,4)
% imagesc(twin1:twin2-1,1:num_units,all_mean_firing_rates2(ii1,:));
% title('Sorted 2nd half-Smoothed By 1st Half')
% colormap('jet')
% axis square
% colorbar
% caxis([pct25 pct975])

corrpeak_95 = prctile(corr_peaks,95);
corr50_95 = prctile(corr50(:),95);
corr_index_95 = prctile(corr_index,95);

figure
subplot(1,2,1)
hist(corr50(:),50)
box off
xlabel('Corr_{1/2}')
ylabel('# Units')
title(['Medan Corr_{1/2} = ' num2str(mean(corr50(:)),3) ', Corr95% = ' ... 
    num2str(corr50_95,3)])

subplot(1,2,2)
hist(corr_index(:),50)
box off
xlabel('Corr Index{1/2}')
ylabel('# Units')
title(['Medan Corr Index {1/2} = ' num2str(mean(corr50(:)),3) ', Corr Index 95% = ' ... 
    num2str(corr_index_95,3)])
