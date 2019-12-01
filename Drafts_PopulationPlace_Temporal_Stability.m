
clar %clear,clc

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size


%---Temporal Stability/Cross Validation
num_shuffs = 1000;
all_out2in = cell(1,100);
all_out2in_nonplace = cell(1,250);
n_out2in = NaN(1,100);
n_out2in_nonplace = NaN(1,250);

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = cell(1,100); %place cell unit names
all_non_place_cell_unit_names = cell(1,250);

monkeys = {'Vivian','Tobii'};
figure_dir = {};
place_cell_ind = 1;
nonplace_cell_ind = 1;
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = 'P:\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalysis\PW Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalysis\PW Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = 'P:\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalysis\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
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
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isnan(stats_across_tasks(1,unit)) %no peak time
                continue
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                if sum(in_out{unit} <= 2) < 0%200
                    continue
                end
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names{place_cell_ind} = [task_file(1:end-11) '_' unit_stats{1,unit}];
                n_out2in(place_cell_ind) = sum(in_out{unit} <= 2);%number of fixaitons out->in
                
                %firing rate out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} <= 2,:); %get spike trains
                [~,all_out2in{place_cell_ind}] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                place_cell_ind = place_cell_ind+1;
            else
                
                if sum(in_out{unit} <= 2) < 0
                    continue
                elseif isnan(stats_across_tasks(1,unit)) %no peak time
                    continue
                end
                
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_nonplace_cell_unit_names{nonplace_cell_ind} = [task_file(1:end-11) '_' unit_stats{1,unit}];
                n_out2in_nonplace(nonplace_cell_ind) = sum(in_out{unit} <= 2);%number of fixaitons out->in
                
                %firing rate out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} <= 2,:); %get spike trains
                [~,all_out2in_nonplace{nonplace_cell_ind}] = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                nonplace_cell_ind = nonplace_cell_ind+1;
            end
        end
    end
end
place_cell_ind = place_cell_ind-1;
nonplace_cell_ind = nonplace_cell_ind-1;
%% For Place Cells

observed_corr50 = NaN(1,place_cell_ind);
shuff_corr50 = NaN(num_shuffs,place_cell_ind);

all_shuff_firing_rates = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates{i,shuff} = NaN(place_cell_ind,twin1+twin2);
    end
end

for unit = 1:place_cell_ind
    num_trials = size(all_out2in{unit},1);
    half_trials = floor(num_trials/2);
    
    
    %     Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
    %     Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
    Dmean1 = mean(all_out2in{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in{unit}(2:2:end,:));
    observed_corr50(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        shuff_firing = circshift_row(all_out2in{unit});
        
        %     Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
        %     Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
        Dmean1 = mean(shuff_firing(1:2:end,:));
        Dmean2 = mean(shuff_firing(2:2:end,:));
        
        all_shuff_firing_rates{1,shuff}(unit,:) = Dmean1;
        all_shuff_firing_rates{2,shuff}(unit,:) = Dmean2;
        
        shuff_corr50(shuff,unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    end
end
%% For Non-place cells

observed_corr50_nonplace = NaN(1,place_cell_ind);
shuff_corr50_nonplace = NaN(num_shuffs,place_cell_ind);

all_shuff_firing_rates_nonplace = cell(2,num_shuffs);
for i = 1:2
    for shuff = 1:num_shuffs
        all_shuff_firing_rates_nonplace{i,shuff} = NaN(nonplace_cell_ind,twin1+twin2);
    end
end


for unit = 1:nonplace_cell_ind
    num_trials = size(all_out2in_nonplace{unit},1);
    half_trials = floor(num_trials/2);
    
    
    %     Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
    %     Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
    Dmean1 = mean(all_out2in_nonplace{unit}(1:2:end,:));
    Dmean2 = mean(all_out2in_nonplace{unit}(2:2:end,:));
    observed_corr50_nonplace(unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    
    for shuff = 1:num_shuffs
        
        shuff_firing = circshift_row(all_out2in_nonplace{unit});
        
        
        %     Dmean1 = mean(all_firing_rates{unit}(1:half_trials,:));
        %     Dmean2 = mean(all_firing_rates{unit}(half_trials+1:end,:));
        Dmean1 = mean(shuff_firing(1:2:end,:));
        Dmean2 = mean(shuff_firing(2:2:end,:));
        
        all_shuff_firing_rates_nonplace{1,shuff}(unit,:) = Dmean1;
        all_shuff_firing_rates_nonplace{2,shuff}(unit,:) = Dmean2;
        
        shuff_corr50_nonplace(shuff,unit) = corr(Dmean1',Dmean2','row','pairwise','type','Spearman');
    end
end
%% Get % of Units with Significant Temporal Stability
prctile_corr50 = NaN(1,place_cell_ind);
for unit = 1:place_cell_ind
    prctile_corr50(unit) = 100*sum(observed_corr50(unit) > shuff_corr50(:,unit))/num_shuffs;
end

prctile_corr50_non_place = NaN(1,nonplace_cell_ind);
for unit = 1:nonplace_cell_ind
    prctile_corr50_non_place(unit) = 100*sum(observed_corr50_nonplace(unit) > shuff_corr50_nonplace(:,unit))/num_shuffs;
end


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
    [~,i1] = sort(mx1);
    
    [~,mx2] = max(firining_rates2,[],2);
    [~,i2] = sort(mx2);
    
    shuff_peak_corrs(shuff) = corr(mx1,mx2,'row','pairwise','type','Spearman');
    shuff_order_corr(shuff) = corr(i1,i2,'row','pairwise','type','Spearman');
    shuff_order_distance(shuff) = mean(abs(mx1-mx2));
end
%% Get Observed of  Place Cell Orders Correlation
firing_rates1 = NaN(place_cell_ind,twin1+twin2);
firing_rates2 = NaN(place_cell_ind,twin1+twin2);
for unit = 1:place_cell_ind
    num_trials = size(all_out2in{unit},1);
    half_trials = floor(num_trials/2);
    
    Dmean1 = mean(all_out2in{unit}(1:half_trials,:));
    Dmean2 = mean(all_out2in{unit}(half_trials+1:end,:));
    
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

%%
figure
subplot(2,2,1)
imagesc(firing_rates1(i1,:))
caxis([-0.2 1])
title('Odd Sorted')

subplot(2,2,2)
imagesc(firing_rates2(i2,:))
caxis([-0.2 1])
title('Even Sorted')

subplot(2,2,3)
imagesc(firing_rates1(i2,:))
caxis([-0.2 1])
title('Odd Sorted by Even')

subplot(2,2,4)
imagesc(firing_rates2(i1,:))
caxis([-0.2 1])
title('Even Sorted by Odd')

%% Get Distribution of Shuffled Non-Place Cell Orders Correlations

shuff_order_corr_nonplace = NaN(1,num_shuffs);
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
    [~,i1] = sort(mx1);
    
    [~,mx2] = max(firining_rates2,[],2);
    [~,i2] = sort(mx2);
    
    shuff_order_corr_nonplace(shuff) = corr(i1,i2,'row','pairwise','type','Spearman');
end
%% Get Observed of  Non-Place Cell Orders Correlation
firing_rates1 = NaN(nonplace_cell_ind,twin1+twin2);
firing_rates2 = NaN(nonplace_cell_ind,twin1+twin2);
for unit = 1:nonplace_cell_ind
    num_trials = size(all_out2in_nonplace{unit},1);
    half_trials = floor(num_trials/2);
    
    Dmean1 = mean(all_out2in_nonplace{unit}(1:half_trials,:));
    Dmean2 = mean(all_out2in_nonplace{unit}(half_trials+1:end,:));
    
    Dmean1 = Dmean1-mean(Dmean1(1:twin1));
    Dmean1 = Dmean1/max(Dmean1);
    
    Dmean2 = Dmean2-mean(Dmean2(1:twin1));
    Dmean2 = Dmean2/max(Dmean2);
    
    firing_rates1(unit,:) = Dmean1;
    firing_rates2(unit,:) = Dmean2;
end
% firing_rates1 = laundry(firing_rates1);
% firing_rates2 = laundry(firing_rates2);
% if size(firing_rates1,1) ~= size(firing_rates2,1)
%    error('What?')
% end
firing_rates1(isinf(firing_rates1)) = 0;
firing_rates2(isinf(firing_rates2)) = 0;

[~,mx1] = max(firing_rates1,[],2);
[~,i1] = sort(mx1);

[~,mx2] = max(firing_rates2,[],2);
[~,i2] = sort(mx2);

observed_order_corr_nonplace = corr(i1,i2,'row','pairwise','type','Spearman');

%% Print Summary of Results
median_trial_count = 220;%real medain is 203.5,median(n_out2in);
disp('-------------------------------------------------------------------------')
disp(['# Place cells with place field temporal stability (>95%): ' num2str(sum(prctile_corr50 > 95)) ...
    '/' num2str(place_cell_ind)]);
disp(['# Place cells with place field marginal temporal stability (>90%): ' num2str(sum(prctile_corr50 > 90)) ...
    '/' num2str(place_cell_ind)]);
disp(['Median temporal stability: ' num2str(median(observed_corr50),3) ',' ...
    'Mean temporal stability: ' num2str(mean(observed_corr50),3)]);
disp('-------------------------------------------------------------------------')
disp(['# Place Cells with temporal stability (> 95%) with more data: ' ...
    num2str(sum(prctile_corr50 > 95 & n_out2in(1:place_cell_ind) > median_trial_count)) '/'...
    num2str(sum(n_out2in > median_trial_count))])
disp(['# Place Cells with temporal stability (> 95%) with less data: ' ...
    num2str(sum(prctile_corr50 > 95 & n_out2in(1:place_cell_ind) < median_trial_count)) '/'...
    num2str(sum(n_out2in < median_trial_count))])
disp(['Median Corr50 for Place Cells with with more data: ' ...
    num2str(mean(observed_corr50(n_out2in(1:place_cell_ind) > median_trial_count)),3)])
disp(['Median Corr50 for Place Cells with with less data: ' ...
    num2str(mean(observed_corr50(n_out2in(1:place_cell_ind) < median_trial_count)),3)])
disp('-------------------------------------------------------------------------')

disp(['Population Place Temporal stability: ' num2str(observed_order_corr,3) ', greater than '...
    num2str(sum(observed_order_corr > shuff_order_corr)/num_shuffs*100,3) '%'])

disp(['Population Non-Place Temporal stability: ' num2str(observed_order_corr_nonplace,3) ', greater than '...
    num2str(sum(observed_order_corr_nonplace > shuff_order_corr_nonplace)/num_shuffs*100,3) '%'])

disp('-------------------------------------------------------------------------')

disp(['# Non Place cells with place field temporal stability (>95%): ' num2str(sum(prctile_corr50_non_place > 95)) ...
    '/' num2str(nonplace_cell_ind)]);
disp(['# Non Place cells with place field marginal temporal stability (>90%): ' num2str(sum(prctile_corr50_non_place > 90)) ...
    '/' num2str(nonplace_cell_ind)]);
disp(['Median temporal stability: ' num2str(nanmedian(observed_corr50_nonplace),3) ',' ...
    'Mean temporal stability: ' num2str(nanmean(observed_corr50_nonplace),3)]);
disp('-------------------------------------------------------------------------')
disp(['Shuffled Order Distance percentile: ' num2str(100*sum(observed_order_distance < shuff_order_distance)/num_shuffs,3)])