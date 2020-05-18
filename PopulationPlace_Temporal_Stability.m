%Code cleaned re-written from Draft_populationnonplace_Temporal_Stability SDK 5/17/2020.
%Code asks if population-level activity is more structured in time than
%expected by random chance. Code looks at distance between peak-times and
%correlation of peak for even/odd trials.

clar %clear,clc

%---Figure Options---%
plot_individual_cells = true;
set(0,'DefaultFigureVisible','ON');
figure_dir2 = 'C:\Users\seth.koenig\Desktop\New folder\';

%---Import Options---%
num_shuffs = 1000;
task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800;
imageY = 600;


%---Fixation/Saccade Duration parameters---%
min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

only_process_active_units = true; %units with criterion below
t_start = 50; %time before fixation to include in information score calculation, this period is not warped
min_saccades_with_spikes = 0.05;% 5% to remove neurons that have basically no activity since they pass too :(...
%essentially 0.4 Hz threshold if window is 50 ms in wide
min_num_fixations = 100;

%---Unit Information (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = cell(1,100); %place cell unit names
all_nonplace_cell_unit_names = cell(1,250);


%---Load In Place Cell Activity---%
place_cell_ind = 1;
nonplace_cell_ind = 1;
skipped_units = 0;
monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = 'P:\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\sethk\OneDrive\Documents\MATLAB\ViewCellPaperAnalysis\PW Recording Files\';
        figure_dir = 'C:\Users\sethk\OneDrive\Documents\MATLAB\ViewCellPaperAnalysis\PW Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = 'P:\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\sethk\OneDrive\Documents\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
        figure_dir = 'C:\Users\sethk\OneDrive\Documents\MATLAB\ViewCellPaperAnalysis\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,...
            sorting_quality,waveform_count,lfp_quality,comments] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat']);
        
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
        
        %load spatial analysis data
        disp(task_file(1:8))
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spatial_info','spike_times','eyepos','binsize','filter_width')
        
        
        %load Place Cell Fixation Analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{sess},task);
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        for unit =1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                if isnan(stats_across_tasks(1,unit))%no peak detected
                    skipped_units = skipped_units+1;
                    continue
                end
                
                
                %--Determine if should NOT process Data do to sparse sampling---%
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains for out-> in
                spike_counts_out2in = nansum(firing_rate(:,twin1-t_start+1:end));
                if sum(spike_counts_out2in > 0) < min_saccades_with_spikes*size(firing_rate,1) || ...
                        min_num_fixations > size(firing_rate,1)
                    if only_process_active_units
                        skipped_units = skipped_units+1;
                        continue
                    end
                end
                
                %---Store unit name---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names{place_cell_ind} = [task_file(1:end-11) '_' unit_stats{1,unit}];
                
                %---store already smoothed firing rate across units---%
                %smoothing nows save processing time later
                [~,all_out2in{place_cell_ind}] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                place_cell_ind = place_cell_ind+1;
            else
                
                %unlikley to be peaky so ignoring this to get more units,
                %otherwise only have 57
                %                 if isnan(stats_across_tasks(1,unit))%no peak detected
                %                     skipped_units = skipped_units+1;
                %                     continue
                %                 end
                
                %--Determine if should NOT process Data do to sparse sampling---%
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains for out-> in
                spike_counts_out2in = nansum(firing_rate(:,twin1-t_start+1:end));
                if sum(spike_counts_out2in > 0) < min_saccades_with_spikes*size(firing_rate,1) || ...
                        min_num_fixations > size(firing_rate,1)
                    if only_process_active_units
                    skipped_units = skipped_units+1;
                    continue
                    end
                end
                
                %---Store unit name---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_nonplace_cell_unit_names{nonplace_cell_ind} = [task_file(1:end-11) '_' unit_stats{1,unit}];
                
                %---store already smoothed firing rate across units---%
                %smoothing nows save processing time later
                [~,all_out2in_nonplace{nonplace_cell_ind}] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                nonplace_cell_ind = nonplace_cell_ind+1;
            end
        end
    end
end


%do some cleanup
all_out2in_nonplace(nonplace_cell_ind:end) = [];
all_out2in(place_cell_ind:end) = [];

place_cell_ind = place_cell_ind-1;
nonplace_cell_ind = nonplace_cell_ind-1;

%% Individual Unit Stats

%store shuffled firing rates so don't have to do this more than once
all_shuff_firing_rates = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates{i,shuff} = NaN(place_cell_ind,twin1+twin2);
    end
end

%---Get Temporal Stability Measurements For Each Place Cell---%
observed_distance50 = NaN(1,nonplace_cell_ind); %distance between peaks
observed_corr50 = NaN(1,nonplace_cell_ind); %correlation between firing rate curves
shuff_distance50 = NaN(num_shuffs,nonplace_cell_ind);
shuff_corr50 = NaN(num_shuffs,nonplace_cell_ind);

for unit = 1:place_cell_ind
    
    Dmean1 = mean(all_out2in{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in{unit}(2:2:end,:));
    
    %get peak times, not using other more complicated methods
    [~,mx1] = max(Dmean1,[],2);
    [~,mx2] = max(Dmean2,[],2);
    
    observed_distance50(unit) = abs(mx1-mx2); 
    observed_corr50(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        shuff_firing = circshift_row(all_out2in{unit});
        
        Dmean1 = mean(shuff_firing(1:2:end,:));
        Dmean2 = mean(shuff_firing(2:2:end,:));
        
        [~,mx1] = max(Dmean1,[],2);
        [~,mx2] = max(Dmean2,[],2);

        shuff_distance50(shuff,unit) = abs(mx1-mx2); 
        shuff_corr50(shuff,unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
        
        all_shuff_firing_rates{1,shuff}(unit,:) = Dmean1;
        all_shuff_firing_rates{2,shuff}(unit,:) = Dmean2;
    end
end

%get proportion of significant units
place_prctile_distance50 = NaN(1,place_cell_ind);
place_prctile_corr50 = NaN(1,place_cell_ind);
for unit = 1:place_cell_ind
    place_prctile_corr50(unit) = 100*sum(observed_corr50(unit) > shuff_corr50(:,unit))/num_shuffs;
    place_prctile_distance50(unit) = 100*sum(observed_distance50(unit) < shuff_distance50(:,unit))/num_shuffs; 
end

%%
%---Get Temporal Stability Measurements For Each Place Cell---%

%store shuffled firing rates so don't have to do this more than once
all_shuff_firing_rates_nonplace = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates_nonplace{i,shuff} = NaN(nonplace_cell_ind,twin1+twin2);
    end
end


nonplace_observed_distance50 = NaN(1,nonplace_cell_ind); %distance between peaks
nonplace_observed_corr50 = NaN(1,nonplace_cell_ind); %correlation between firing rate curves
nonplace_shuff_distance50 = NaN(num_shuffs,nonplace_cell_ind);
nonplace_shuff_corr50 = NaN(num_shuffs,nonplace_cell_ind);

for unit = 1:nonplace_cell_ind
    
    Dmean1 = mean(all_out2in_nonplace{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in_nonplace{unit}(2:2:end,:));
    
    %get peak times, not using other more complicated methods
    [~,mx1] = max(Dmean1,[],2);
    [~,mx2] = max(Dmean2,[],2);
    
    nonplace_observed_distance50(unit) = abs(mx1-mx2); 
    nonplace_observed_corr50(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        shuff_firing = circshift_row(all_out2in_nonplace{unit});
        
        Dmean1 = mean(shuff_firing(1:2:end,:));
        Dmean2 = mean(shuff_firing(2:2:end,:));
        
        [~,mx1] = max(Dmean1,[],2);
        [~,mx2] = max(Dmean2,[],2);

        nonplace_shuff_distance50(shuff,unit) = abs(mx1-mx2); 
        nonplace_shuff_corr50(shuff,unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
        
        all_shuff_firing_rates_nonplace{1,shuff}(unit,:) = Dmean1;
        all_shuff_firing_rates_nonplace{2,shuff}(unit,:) = Dmean2;
    end
end

%get proportion of significant units
nonplace_prctile_distance50 = NaN(1,nonplace_cell_ind);
nonplace_prctile_corr50 = NaN(1,nonplace_cell_ind);
for unit = 1:nonplace_cell_ind
    nonplace_prctile_corr50(unit) = 100*sum(nonplace_observed_corr50(unit) > nonplace_shuff_corr50(:,unit))/num_shuffs;
    nonplace_prctile_distance50(unit) = 100*sum(nonplace_shuff_distance50(unit) < nonplace_shuff_distance50(:,unit))/num_shuffs; 
end


%significance of plave vs non-place
pvalCorr50 = chiSquareProportionTest(sum(place_prctile_corr50 > 95),place_cell_ind,...
    sum(nonplace_prctile_corr50 > 95),nonplace_cell_ind);
pvalDistance50 = chiSquareProportionTest(sum(place_prctile_distance50 > 95),place_cell_ind,...
    sum(nonplace_prctile_distance50 > 95),nonplace_cell_ind);

%% Get Observed of  Place Cell Orders Correlation
firing_rates1 = NaN(place_cell_ind,twin1+twin2);
firing_rates2 = NaN(place_cell_ind,twin1+twin2);
for unit = 1:place_cell_ind
    
    Dmean1 = mean(all_out2in{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in{unit}(2:2:end,:));
    
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
observed_order_distance = mean(abs(mx1-mx2));

figure
plot(mx1-twin1,mx2-twin2,'.k')
xlabel('Peak Response on Odd Trials (ms)')
ylabel('Peak Response on Even Trials (ms)')
box off
title(['Place Cell Population-level Peak Correlation (\rho) = ' num2str(observed_peak_corrs,3)])

%% Get Distribution of Shuffled Place Cell Orders Correlations

shuff_order_corr = NaN(1,num_shuffs);
shuff_peak_corrs =  NaN(1,num_shuffs);
shuff_order_distance = NaN(1,num_shuffs);
for shuff = 1:num_shuffs
    
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
    [~,mx2] = max(firining_rates2,[],2);
    
    shuff_peak_corrs(shuff) = corr(mx1,mx2,'row','pairwise','type','Spearman');
    shuff_order_distance(shuff) = mean(abs(mx1-mx2));
end

%% Get Observed of  nonplace Cell Orders Correlation
firing_rates1 = NaN(nonplace_cell_ind,twin1+twin2);
firing_rates2 = NaN(nonplace_cell_ind,twin1+twin2);
for unit = 1:nonplace_cell_ind
    
    Dmean1 = mean(all_out2in_nonplace{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in_nonplace{unit}(2:2:end,:));
    
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

nonplace_observed_peak_corrs = corr(mx1,mx2,'row','pairwise','type','Spearman');
nonplace_observed_order_distance = mean(abs(mx1-mx2));

figure
plot(mx1-twin1,mx2-twin2,'.k')
xlabel('Peak Response on Odd Trials (ms)')
ylabel('Peak Response on Even Trials (ms)')
box off
title(['Non-place Cell Population Peak Correlation (\rho) = ' num2str(nonplace_observed_peak_corrs,3)])

nonplace_shuff_order_corr = NaN(1,num_shuffs);
nonplace_shuff_peak_corrs =  NaN(1,num_shuffs);
nonplace_shuff_order_distance = NaN(1,num_shuffs);
for shuff = 1:num_shuffs
    
    firining_rates1 = all_shuff_firing_rates_nonplace{1,shuff};
    firining_rates2 = all_shuff_firing_rates_nonplace{2,shuff};
    
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
    [~,mx2] = max(firining_rates2,[],2);
    
    nonplace_shuff_peak_corrs(shuff) = corr(mx1,mx2,'row','pairwise','type','Spearman');
    nonplace_shuff_order_distance(shuff) = mean(abs(mx1-mx2));
end

%% Print Results
clc
disp('------------------------------------------------------')
disp('Inidividual Unit Stats:')
disp(['Place Cell Correlation: ' num2str(sum(place_prctile_corr50 > 95)) '(' num2str(100*sum(place_prctile_corr50 > 95)/place_cell_ind,3) '%)'])
disp(['Place Cell Peak Distance: ' num2str(sum(place_prctile_distance50 > 95)) '(' num2str(100*sum(place_prctile_distance50 > 95)/place_cell_ind,3) '%)'])
disp(' ')
disp(['Nonplace Cell Correlation: ' num2str(sum(nonplace_prctile_corr50 > 95)) '(' num2str(100*sum(nonplace_prctile_corr50 > 95)/nonplace_cell_ind,3) '%)'])
disp(['Nonplace Cell Peak Distance: ' num2str(sum(nonplace_prctile_distance50 > 95)) '(' num2str(100*sum(nonplace_prctile_distance50 > 95)/nonplace_cell_ind,3) '%)'])
disp(' ')
disp(['Signficance place vs nonplace Correlation: ' num2str(pvalCorr50)])
disp(['Signficance place vs nonplace Distance: ' num2str(pvalDistance50)])
disp('------------------------------------------------------')
disp('')
disp(['Place Cell Population Order-Peak: ' num2str(100*sum(observed_peak_corrs > shuff_peak_corrs)/num_shuffs) '%'])
disp(['Place Cell Population Order-Distance: ' num2str(100*sum(observed_order_distance < shuff_order_distance)/num_shuffs) '%'])
disp('')
disp(['NonPlace Cell Population Order-Peak: ' num2str(100*sum(nonplace_observed_peak_corrs > nonplace_shuff_peak_corrs)/num_shuffs) '%'])
disp(['NonPlace Cell Population Order-Distance: ' num2str(100*sum(nonplace_observed_order_distance < nonplace_shuff_order_distance)/num_shuffs) '%'])


