% Code below creates population summary for Significnat Saccade Direction Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes MRLs (mean resultant vector length) for place and non-Saccade Direction Cells
% 2) Summarize circular non-uniformity p-value (biased by fixation count and firing rate)
% 3) Copies relevant figures for Saccade Direction Cells to summary directory

%Code rechecked by SDK on 1/5/2017 & then re-written and rechecked again 5/3/2017

clar %clear,clc

%where to store spatial analysis figure copies
summary_directory =  'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PopulationFigures\Saccade Direction\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
saccades_with_spikes = [];
smval_deg = 6;%18 %9 degrees std

num_shuffs = 1000;
all_firing_rates = cell(1,100);
cell_ind = 1;
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data');
        
        %get unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        if num_units == 0
            continue
        end
        
        num_trials = length(cfg.trl); %number of trials
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            continue
        end
        
        disp(task_file(1:8))
        if exist([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'])
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
            load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        twinad1 = 200;%presaccade window
        twinad2 = 400;%pos saccade window
        start_window = twinad1-min_fix_dur;%how early before saccade can you look for direction tuning
        end_window = twinad1+min_fix_dur+44;%how late after saccade can you look for direction tuning
        %minimum fixation duration + median saccade duration of 44 ms, some neurons have
        
        num_units = size(unit_stats,2);
        for unit =1:num_units
            if ~isnan(mrls.all_saccades(unit)) %if unit was processed
                if  mrls.all_saccades_shuffled_prctile(unit) > 95
                    
                    window = all_direction_windows{unit};
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%---Grab Saccade Aligned Activity--%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    sac_aligned = saccade_aligned_firing{unit}; %fixation aligned firing rate
                    sac_dirs = saccade_directions{unit}; %saccade directions organized the same way as sac_algined
                    fix_starts = fixation_starts{unit};
                    fix_ends = fixation_ends{unit};
                    
                    %---Remove Counfounding Eye Movements for Spatially Modulated Neurons---%
                    %take eye data from fixations in2in or out2out since in2out or out2in
                    %could be biased by field location creating artificial direction tuning
                    % but only do this if spatial (both criterion)
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        sac_aligned = sac_aligned(sac_in_out{unit} == 2 | sac_in_out{unit} == 4,:); %fixations in2in or out2out
                        sac_dirs = sac_dirs(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for fixations in2in or out2out only
                        fix_starts = fix_starts((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
                        fix_ends = fix_ends((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
                    end
                    
                    %---Remove Fixations that are too short in duration---%
                    window_width = length(all_direction_windows{unit});
                    window_end = all_direction_windows{unit}(end);
                    window_start = all_direction_windows{unit}(1);
                    if window_end > end_window
                        window_end = all_direction_windows{unit}(end)-twinad1;
                        %remove fixations shorter than window as these could be
                        %contaminated by the next saccade
                        fixations_too_short = find(fix_ends < window_end);
                        sac_dirs(fixations_too_short) = [];
                        sac_aligned(fixations_too_short,:) = [];
                    end
                    if window_start < twinad1
                        fixations_too_short = find((twinad1+fix_starts) > window_start);
                        sac_dirs(fixations_too_short) = [];
                        sac_aligned(fixations_too_short,:) = [];
                    end
                    
                    [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_aligned,window,sac_dirs);
                    binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
                    [estimated_prefered_direction,prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
                        select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,smval_deg,bin_deg);
                    
                  
                    [~,all_firing_rates{cell_ind}] = nandens3(sac_aligned(prefered_dirs,:),smval,1000);

                    cell_ind = cell_ind + 1;
                end
            end
        end
    end
end
%%
all_firing_rates(cell_ind:end) = [];
cell_ind = cell_ind-1;
num_sacs = cellfun('size',all_firing_rates,1);

%% Get Observed  Orders Correlation
firing_rates1 = NaN(cell_ind,twin1+twin2);
firing_rates2 = NaN(cell_ind,twin1+twin2);
for unit = 1:cell_ind
    num_trials = size(all_firing_rates{unit},1);
    half_trials = floor(num_trials/2);
    
    Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
    Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
    
    Dmean1 = Dmean1-mean(Dmean1(1:twin1));
    Dmean1 = Dmean1/max(Dmean1);
    
    Dmean2 = Dmean2-mean(Dmean2(1:twin1));
    Dmean2 = Dmean2/max(Dmean2);
    
    firing_rates1(unit,:) = Dmean1;
    firing_rates2(unit,:) = Dmean2;
end

[~,mx1] = max(firing_rates1,[],2);
[~,i1] = sort(mx1);

[~,mx2] = max(firing_rates2,[],2);
[~,i2] = sort(mx2);

observed_peak_corrs = corr(mx1,mx2,'row','pairwise','type','Spearman');
observed_order_corr = corr(i1,i2,'row','pairwise','type','Spearman');
observed_order_distance = mean(abs(mx1-mx2));

figure
plot(mx1-twin1,mx2-twin1,'.k')
title(['Observed Peak Corrs Even/Odd Trials: \rho = ' num2str(observed_peak_corrs,3)])
xlabel('Time from Saccde Onset (ms)')
ylabel('Time from Saccade Onset (ms)')
axis square
xlim([-twin1 twin2])
ylim([-twin1 twin2])
box off
%%
observed_corr50 = NaN(1,cell_ind);
shuff_corr50 = NaN(num_shuffs,cell_ind);

all_shuff_firing_rates = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates{i,shuff} = NaN(cell_ind,twin1+twin2);
    end
end

for unit = 1:cell_ind
    num_trials = size(all_firing_rates{unit},1);
    half_trials = floor(num_trials/2);
    
%     Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
%     Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
    Dmean1 = mean(all_firing_rates{unit}(1:2:end,:));
    Dmean2 = mean(all_firing_rates{unit}(2:2:end,:));

    observed_corr50(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        shuff_firing = circshift_row(all_firing_rates{unit});
%         
%         Dmean1 = mean(shuff_firing(1:half_trials,:));
%         Dmean2 = mean(shuff_firing(half_trials+1:end,:));
        
        Dmean1 = mean(shuff_firing(1:2:end,:));
        Dmean2 = mean(shuff_firing(2:2:end,:));
        
        all_shuff_firing_rates{1,shuff}(unit,:) = Dmean1;
        all_shuff_firing_rates{2,shuff}(unit,:) = Dmean2;
        
        shuff_corr50(shuff,unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    end
end
%
%% Get Distribution of Shuffled Orders Correlations

shuff_order_corr = NaN(1,num_shuffs);
shuff_peak_corrs =  NaN(1,num_shuffs);
shuff_order_distance = NaN(1,num_shuffs);
for shuff = 1:num_shuffs
    
    %
    %         Dmean1 = mean(shuff_firing(1:half_trials,:));
    %         Dmean2 = mean(shuff_firing(half_trials+1:end,:));
    Dmean1 = mean(shuff_firing(1:2:end,:));
    Dmean2 = mean(shuff_firing(2:2:end,:));
    
    
    firining_rates1 = all_shuff_firing_rates{1,shuff};
    firining_rates2 = all_shuff_firing_rates{2,shuff};
    
    for unit = 1:size(firining_rates1)
        fr = firining_rates1(unit,:);
        fr = fr-mean(fr(1:twin1));
        fr = fr/max(fr);
        firining_rates1(unit,:) = fr;
        
        fr = firining_rates2(unit,:);
        fr = fr-mean(fr(1:twin1));
        fr = fr/max(fr);
        firining_rates2(unit,:) = fr;
    end
    
    
    [~,mx1] = max(firining_rates1,[],2);
    [~,i1] = sort(mx1);
    
    [~,mx2] = max(firining_rates2,[],2);
    [~,i2] = sort(mx2);
    
    shuff_peak_corrs(shuff) = corr(mx1,mx2,'row','pairwise','type','Spearman');
    shuff_order_corr(shuff) = corr(i1,i2,'row','pairwise','type','Spearman');
    shuff_order_distance(shuff) = mean(abs(mx1-mx2));
end

%% Get % of Units with Significant Temporal Stability
prctile_corr50 = NaN(1,cell_ind);
for unit = 1:cell_ind
    prctile_corr50(unit) = 100*sum(observed_corr50(unit) > shuff_corr50(:,unit))/num_shuffs;
end

prctile_peak_corr = 100*sum(observed_peak_corrs > shuff_peak_corrs)/num_shuffs;
prctile_order_corr = 100*sum(observed_order_corr > shuff_order_corr)/num_shuffs; 
prctile_distance = 100*sum(observed_order_distance < shuff_order_distance)/num_shuffs;
%%
%% Print Summary of Results
median_trial_count = median(num_sacs);%200
disp('-------------------------------------------------------------------------')
disp(['# Saccade Direction Cells with temporal stability (>95%): ' num2str(sum(prctile_corr50 > 95)) ...
    '/' num2str(cell_ind)]);
disp(['# Saccade Direction Cells with marginal temporal stability (>90%): ' num2str(sum(prctile_corr50 > 90)) ...
    '/' num2str(cell_ind)]);
disp(['Median temporal stability: ' num2str(median(observed_corr50),3) ',' ...
    'Mean temporal stability: ' num2str(mean(observed_corr50),3)]);
disp('-------------------------------------------------------------------------')
disp(['# Saccade Direction Cells with temporal stability (> 95%) with more data: ' ...
    num2str(sum(prctile_corr50 > 95 & num_sacs(1:cell_ind) > median_trial_count)) '/'...
    num2str(sum(num_sacs > median_trial_count))])
disp(['# Saccade Direction Cells with temporal stability (> 95%) with less data: ' ...
    num2str(sum(prctile_corr50 > 95 & num_sacs(1:cell_ind) < median_trial_count)) '/'...
    num2str(sum(num_sacs < median_trial_count))])
disp(['Median Corr50 for Saccade Direction Cells with with more data: ' ...
    num2str(mean(observed_corr50(num_sacs(1:cell_ind) > median_trial_count)),3)])
disp(['Median Corr50 for Saccade Direction Cells with with less data: ' ...
    num2str(mean(observed_corr50(num_sacs(1:cell_ind) < median_trial_count)),3)])
disp('-------------------------------------------------------------------------')
%%
disp(['Population Direction Temporal stability: ' num2str(observed_order_corr,3) ', greater than '...
    num2str(sum(observed_order_corr > shuff_order_corr)/num_shuffs*100,3) '%'])
