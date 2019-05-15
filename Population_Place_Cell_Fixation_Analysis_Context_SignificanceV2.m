% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes spatial correlations for place cell and non-place cells
% 2) Calculates fixation aligned firing rate curves for place cells
% 3) Calculates eye coverage and place field coverage for place cells
% 4) Tracks AP location, unit counts, and which monkey (not currently used)
% 5) Contextual differences between list and sequence task
% 6) Copies relevant figures for place cells to summary directory

%Streamlined from Population_Place_Cell_Fixation_Analysis_Context_Significance
%updated on 2/12/18 SDK

clar %clear,clc
set(0,'DefaultFigureVisible','OFF');

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Spatial Analysis\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s
place_cell_AP_location = []; %AP location of recorded place cell
place_cell_subregion = []; %e.g. CA3 or CA1
non_place_cell_subregion = [];

n_in = []; %number of fixations ?->in
n_out2in = [];%number of fixaitons out->in

%---Spatial Correlation Values for all cells---%
all_place_cell_spatial_corrs = []; %place cell spatial correlations
all_place_cell_spatial_skagg_percents= []; %place cell skagg information score percentiles
all_non_place_cell_spatial_corrs = [];%non-place cell spatial correlations
all_non_place_cell_skagg_percents = [];%non-place cell skagg information score percentiles
non_place_skagg_names = [];

%---place field properties---%
place_field_area = []; %area of place field
place_coverage = {}; %place field locations
coverage = {}; %eye data coverage for place cells

%---Place Cell Firing Rate Curves for List---%
all_in_rates = [];%fixations out-> in
all_in_minus_out_rates = [];%out-> in minus out-> out
all_out_rates = []; %fixaiton out-> out normalized by max of out->out
all_fixation_rates = [];%all fixaions normalized by max of all fixations
all_list_peak_times = [];%time of peak firing rate of place cells for out-> in fixations
sig_p_list = []; %whether fixation firing rates were reliably different for in vs out for list
all_peak_list = []; %peak firing rate during list
all_list_peak_time = [];%peak firing time list
all_whole_session_firing_rate = [];
putative_excite_inhibit = [];

%---Firing Rate Curve Properties for Sequence Task---%
sig_p_seq = []; %whether firing rates were significantly different for in vs out for seq
all_peak_seq = []; %peak firing rate during seq
all_seq_peak_times = [];%peak firing time seq
all_seq_firing_curves = [];
total_sig_time = [];


%---Firing Rates Between Sequence and List Tasks---%
all_context_gain = []; %change in firing rate (list_fr-seq_fr)/list_fr
all_context_gain2 = [];%gain list_fr/seq_fr
task_context_sig_count = [];%time curves
fr_task_context_count1 = [];%windowed
fr_task_context_count2 = [];%peak rate
fr_context_faster_slower = [];%1 faster, 0 no change, -1, slower
numshuffs = 10000;

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
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
            'stability_attribute','cfg','hdr','data','whole_session_mean_firing_rate','excitatory_inhibitory');
        
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
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                
                %---Misc Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s for monkey
                place_cell_AP_location = [place_cell_AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                all_whole_session_firing_rate = [all_whole_session_firing_rate whole_session_mean_firing_rate(unit)];
                putative_excite_inhibit = [putative_excite_inhibit excitatory_inhibitory(unit)];
                
                subregion = session_data{sess}.subregion;
                vals = textscan(subregion,'%s','Delimiter',',');
                if length(vals{1}) == 1
                    place_cell_subregion = [place_cell_subregion vals{1}];
                else
                    chan = str2double(unit_stats{1,unit}(6));
                    place_cell_subregion = [place_cell_subregion vals{1}(chan)];
                end
                
                %---Spatial Correlation Values for Place cells---%
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_place_cell_spatial_skagg_percents = [all_place_cell_spatial_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                
                
                %---place field properties---%
                H = define_spatial_filter(filter_width); %spatial filter
                trial_data{1} = eyepos{unit}; %eye data
                trial_data{2} = spike_times{unit}; %spike times
                [~,timemaps] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
                coverage = [coverage {timemaps}];  %eye data coverage for place cells
                place_coverage = [place_coverage all_place_field_matrix(unit)]; %place field location
                place_field_area = [place_field_area area(unit)]; %size of place field relative to image size
                
                
                %---Place Cell Firing Rate Curves for List---%
                n_in = [n_in sum(in_out{unit} == 1 | in_out{unit} == 2)]; %number of fixations ?->in
                n_out2in = [n_out2in sum(in_out{unit} == 1)];%number of fixaitons out->in
                
                
                %firing rate out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains
                in_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                in_curve2 = in_curve;
                in_curve = in_curve-nanmean(in_curve(1:twin1)); %remove base line
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    in_curve = in_curve/in_curve(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    in_curve = in_curve/max(in_curve);
                end
                all_in_rates = [all_in_rates; in_curve];%fixations out-> in
                all_peak_list = [all_peak_list stats_across_tasks(2,unit)]; %peak firing rate during list
                all_list_peak_times = [all_list_peak_times stats_across_tasks(1,unit)];%time of peak firing rate of place cells for out-> in fixations
                
                %firing rate out-> out
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %get spike trains
                out_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                out_curve2 = out_curve;
                out_curve = out_curve-nanmean(out_curve(1:twin1)); %remove base line
                out_curve = out_curve/nanmax(out_curve); %not sure what peak would be so normalize by max
                all_out_rates = [all_out_rates; out_curve];%fixation out->out
                
                in_minus_out = in_curve2-out_curve2;
                %in_minus_out = in_minus_out-mean(in_minus_out(1:twin1));
                in_minus_out = in_minus_out/max(in_minus_out);
                all_in_minus_out_rates = [all_in_minus_out_rates; in_minus_out];
                
                %firing rate for all fixations
                firing_rate = list_fixation_locked_firing{unit};%get spike trains
                firing_rate = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                firing_rate = firing_rate-nanmean(firing_rate(1:twin1)); %remove base line
                firing_rate = firing_rate/nanmax(firing_rate);%not sure what peak would be so normalize by max
                all_fixation_rates = [all_fixation_rates; firing_rate]; %all fixations
                
                
                %find when firing rate for in field fixations is higher than out of field
                %for fixation out->in vs out->out
                if sum(list_sig_times{1,unit}) > 0
                    pass = 1;
                else
                    pass = 0;
                end
                
                %find when firing rate for in field fixations is higher than out of field
                %for fixation ?->in vs ?->out
                if sum(list_sig_times{2,unit}) > 0
                    pass2 = 1;
                else
                    pass2 = 0;
                end
                
                if pass == 1 && pass2 == 1 %both show sig difference
                    sig_p_list = [sig_p_list 3];
                elseif pass == 1 %just out->in vs out-> out shows difference
                    sig_p_list = [sig_p_list 1];
                elseif pass2 == 1 %just ?->in vs ?-> out shows difference
                    sig_p_list = [sig_p_list 2];
                else %no difference shown :(
                    sig_p_list = [sig_p_list 0];
                end
                
                %---Firing Rate Curve Properties for Sequence Task---%
                if ~isnan(stats_across_tasks(4,unit))
                    all_peak_seq = [all_peak_seq stats_across_tasks(4,unit)]; %peak firing rate during seq
                end
                all_seq_peak_times = [all_seq_peak_times stats_across_tasks(3,unit)];%peak firing time seq
                
                
                %---Firing Rates Between Sequence and List Tasks---%
                %change in firing rate (list_fr-seq_fr)/list_fr
                all_context_gain = [all_context_gain (stats_across_tasks(2,unit)-stats_across_tasks(4,unit))/stats_across_tasks(2,unit)];
                all_context_gain2 = [all_context_gain2 stats_across_tasks(2,unit)/stats_across_tasks(4,unit)];%gain list_fr/seq_fr
                
                %indexed weirdly so have to do it this way
                all_which_sequence = laundry(all_which_sequence);
                fixation_firing = [];
                for c = 1:4
                    for seq = 1:2
                        fixation_firing = [fixation_firing; sequence_fixation_locked_firing{c,unit}(all_which_sequence{unit}(trial_nums{c,unit}) == seq,:)];
                    end
                end
                
                new_seq_sig = [];
                if ~isempty(in_out_sequence{unit}) %at least 1 item in field and 1 item out of field
                    if  isempty(sequence_sig_times{unit}) %no significant difference
                        sig_p_seq = [sig_p_seq 0];
                        all_seq_firing_curves = [all_seq_firing_curves; NaN(1,twin1+twin2)];
                    else
                        seq_in_out = in_out_sequence{unit};
                        in_curve = nandens(fixation_firing(seq_in_out == 1,:),smval,'gauss',Fs,'nanflt');%smoothed firing rates for in field fixations
                        out_curve = nandens(fixation_firing(seq_in_out == 0,:),smval,'gauss',Fs,'nanflt');%smoothed firing rates for out of field fixations
                        pos_ind = (in_curve-out_curve) > 0; %find expected in field > out of field
                        pos_ind(sequence_sig_times{unit} == 0) = 0; %remove times that are not significant
                        pos_ind = find(pos_ind);
                        neg_ind =  (in_curve-out_curve) < 0 ; %find unexpected in field < out of field
                        neg_ind(sequence_sig_times{unit} == 0) = 0; %remove times that are not significant
                        neg_ind = find(neg_ind);
                        
                        if isnan(stats_across_tasks(4,unit))
                            all_peak_seq = [all_peak_seq max(in_curve)]; %peak firing rate during seq
                        end
                        
                        if mean(in_curve(1:twin1)) > mean(in_curve(twin1:end))
                            in_curve = in_curve-nanmean(in_curve(1:twin1));%maintains that first part is higher
                        else
                            in_curve = in_curve-nanmean(in_curve(twin1+1:end));
                        end
                        in_curve = in_curve/max(in_curve);
                
%                         in_curve = in_curve-nanmean(in_curve(1:twin1));
%                         in_curve = in_curve/max(in_curve);
                        
                        all_seq_firing_curves = [all_seq_firing_curves; in_curve];
                        
                        %find contiguous positive significant regions that
                        %at least 2std in length
                        rmv = [];
                        gaps_pos = findgaps(pos_ind); %find breaks
                        total_pos = 0;%total postive time
                        if ~isempty(gaps_pos)
                            for g = 1:size(gaps_pos,1)
                                gp = gaps_pos(g,:);
                                gp(gp == 0) = [];
                                if length(gp) < 50;%1.5*smval %3  standard deviations
                                    rmv = [rmv g];
                                else
                                    total_pos = total_pos+length(gp);
                                end
                            end
                        end
                        if ~isempty(rmv)
                            gaps_pos(rmv,:) = []; %remove time points that are too short
                        end
                        
                        %find contiguous positive significant regions that
                        %at least 2std in length
                        rmv = [];
                        gaps_neg = findgaps(neg_ind);  %find breaks
                        total_neg = 0;%total negative time
                        if ~isempty(gaps_neg)
                            for g = 1:size(gaps_neg,1)
                                gp = gaps_neg(g,:);
                                gp(gp == 0) = [];
                                if length(gp) < 50;%1.5*smval %3  standard deviations
                                    rmv = [rmv g]; %remove time points that are too short
                                else
                                    total_neg = total_neg+length(gp);
                                end
                            end
                        end
                        if ~isempty(rmv)
                            gaps_neg(rmv,:) = [];
                        end
                        
                        if total_pos > 0 && total_neg == 0 %if only expected result
                            sig_p_seq = [sig_p_seq 1];
                            new_seq_sig =gaps_pos(1:end);
                            new_seq_sig(new_seq_sig == 0) = [];
                        elseif total_pos == 0 && total_neg == 0 %no result
                            sig_p_seq = [sig_p_seq 0];
                            new_seq_sig = [];
                        elseif total_pos > 0 && total_neg > 0 %need to probe further since mixed result
                            
                            %find if negative result is before the start of
                            %the fixation since this could be from ITI or
                            %from another fixation (which could be
                            %influenced by this one). Both cases have been
                            %visually confirmed!
                            rmv = [];
                            for g = 1:size(gaps_neg,1)
                                gp = gaps_neg(g,:);
                                gp(gp == 0) = [];
                                if any(gp < twin1)  %contigous negative started before fixation
                                    rmv = [rmv g];
                                end
                            end
                            gaps_neg(rmv,:) = [];
                            total_neg = sum(sum(gaps_neg > 0));
                            
                            if isempty(gaps_neg) %so only negative time is before fixation so call postive effect
                                sig_p_seq = [sig_p_seq 1];
                                new_seq_sig =gaps_pos(1:end);
                                new_seq_sig(new_seq_sig == 0) = [];
                            elseif total_pos > 2*total_neg %much more positive than neg so call positive effect
                                sig_p_seq = [sig_p_seq 1];
                                new_seq_sig = gaps_pos(1:end);
                                new_seq_sig(new_seq_sig == 0) = [];
                                disp('much more positive')
                            elseif 2*total_pos < total_neg %much more negative than positive so call negative effect
                                sig_p_seq = [sig_p_seq -1];
                                new_seq_sig = gaps_neg(1:end);
                                new_seq_sig(new_seq_sig == 0) = [];
                                disp('much more negative')
                            else %not sure so remove
                                sig_p_seq = [sig_p_seq NaN];
                                disp('unsure so removed')
                                new_seq_sig = [];
                            end
                        elseif total_neg > 0 %if only negative result
                            sig_p_seq = [sig_p_seq -1];
                            new_seq_sig = gaps_neg(1:end);
                            new_seq_sig(new_seq_sig == 0) = [];
                        else
                            error('What else could happen')
                        end
                    end
                else
                    sig_p_seq = [sig_p_seq NaN]; %no shapes in and out of field so keep parallel structure
                    all_peak_seq = [all_peak_seq NaN];
                    all_seq_firing_curves = [all_seq_firing_curves; NaN(1,twin1+twin2)];
                end
                
                if~isnan(sig_p_seq(end))
                    if isnan(all_context_gain2(end)) %can happen if no definable peaks so set as ratio of maximum firing rates
                        all_context_gain2(end) =  all_peak_list(end)/all_peak_seq(end);
                    elseif all_context_gain2(end) ~= all_peak_list(end)/all_peak_seq(end);
                        error('Contextual gain 2 should be list peak/sequence peak')
                    end
                    
                    if isnan(all_context_gain(end)) %can happen if no definable peaks so set as ratio of maximum firing rates
                        all_context_gain(end) = (all_peak_list(end)-all_peak_seq(end))/all_peak_list(end);
                    elseif  all_context_gain(end) ~= (all_peak_list(end)-all_peak_seq(end))/all_peak_list(end)
                        error('Context gain should be the change in firing rate')
                    end
                    
                    %indexed weirdly so have to do it this way
                    fixation_firing = [];
                    for c = 1:4
                        for seq = 1:2
                            fixation_firing = [fixation_firing; sequence_fixation_locked_firing{c,unit}(all_which_sequence{unit}(trial_nums{c,unit}) == seq,:)];
                        end
                    end
                    sequence_firing_rate = fixation_firing(find(in_out_sequence{unit} == 1),:);
                    out_sequence_firing_rate = fixation_firing(find(in_out_sequence{unit} == 0),:);
                    list_in_firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains
                    list_out_firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %get spike trains
                    
                    %% in field across tasks
                    [~,list_firing_rate_curves] = nandens(list_in_firing_rate,smval,'gauss',Fs,'nanflt'); %trial by trial smoothing
                    [~,seq_firing_rate_curves] = nandens(sequence_firing_rate,smval,'gauss',Fs,'nanflt'); %trial by trial smoothing
                    both_task_curves = [list_firing_rate_curves; seq_firing_rate_curves];
                    both_task_ind = [ones(size(list_firing_rate_curves,1),1); zeros(size(seq_firing_rate_curves,1),1)];
                    all_curves = NaN(numshuffs,twin1+twin2);
                    parfor shuff = 1:numshuffs;
                        ind = randperm(length(both_task_ind));
                        shuff_task = both_task_ind(ind);
                        shuff_task1_curve = nanmean(both_task_curves(shuff_task == 1,:));
                        shuff_task2_curve = nanmean(both_task_curves(shuff_task == 0,:));
                        all_curves(shuff,:) = shuff_task1_curve-shuff_task2_curve;
                    end
                    observed_diff = nanmean(list_firing_rate_curves)-nanmean(seq_firing_rate_curves); %observed difference in firing rate
                    [~,task_sig_times] = cluster_level_statistic(observed_diff,all_curves,2,smval); %multiple comparision corrected significant indeces
                    
                    %%
                    if ~isempty(new_seq_sig)
                        nss = zeros(1,twin1+twin2);
                        nss(new_seq_sig) = 2;
                    else
                        nss = list_sig_times{1,unit};
                    end
                    
                    sig = 0;
                    gaps = findgaps(find(task_sig_times));
                    good_gp = [];
                    for g = 1:size(gaps,1)
                        gp = gaps(g,:);
                        gp(gp < twin1/2) = [];
                        gp(gp == 0) = [];
                        if length(gp) > 50
                            sig = 1;
                            good_gp = [good_gp gp];
                        end
                    end
                    task_sig_times = zeros(1,twin1+twin2);
                    task_sig_times(good_gp) = 1;
                    
                    if sig == 1
                        task_context_sig_count = [task_context_sig_count 1];
                    else
                        task_context_sig_count =  [task_context_sig_count 0];
                    end
                    
                    
                    %look for differences in firing rate during identified window
                    if sum(list_sig_times{1,unit}) == 0;%happens very rarely
                        window_list = nanmean(list_firing_rate_curves(:,list_sig_times{2,unit} == 2)');%use all in vs all out
                    else
                        window_list = nanmean(list_firing_rate_curves(:,list_sig_times{1,unit} == 2)');%use out2in vs out2out
                    end
                    window_seq = nanmean(seq_firing_rate_curves(:,nss == 2)');
                    
                    
                    %                     window_diff = nanmean(window_list)-nanmean(window_seq);
                    %                     window_all_values = [window_list window_seq];
                    %                     task_id = [ones(1,length(window_list)) 2*ones(1,length(window_seq))];
                    %                     shuff_window_diff = NaN(1,numshuffs);
                    %                     parfor shuff = 1:numshuffs;
                    %                         ind = randperm(length(task_id));
                    %                         shuff_task = task_id(ind);
                    %                         shuff_task1_curve = nanmean(window_all_values(shuff_task == 1));
                    %                         shuff_task2_curve = nanmean(window_all_values(shuff_task == 2));
                    %                         shuff_window_diff(shuff) = shuff_task1_curve-shuff_task2_curve;
                    %                     end
                    %
                    %                     %look for peak firing rate differenes
                    %                     max_list = nanmax(list_firing_rate_curves');
                    %                     max_seq = nanmax(seq_firing_rate_curves');
                    %                     max_diff = nanmean(max_list)-nanmean(max_seq);
                    %                     max_all_values = [max_list max_seq];
                    %                     task_id = [ones(1,length(max_list)) 2*ones(1,length(max_seq))];
                    %                     shuff_max_diff = NaN(1,numshuffs);
                    %                     parfor shuff = 1:numshuffs;
                    %                         ind = randperm(length(task_id));
                    %                         shuff_task = task_id(ind);
                    %                         shuff_task1_curve = nanmean(max_all_values(shuff_task == 1));
                    %                         shuff_task2_curve = nanmean(max_all_values(shuff_task == 2));
                    %                         shuff_max_diff(shuff) = shuff_task1_curve-shuff_task2_curve;
                    %                     end
                    %%
                    %                     t = -twin1:twin2-1;
                    %                     figure
                    %                     subplot(2,3,1)
                    %                     dofill(t,list_in_firing_rate,'red',1,smval);
                    %                     hold on
                    %                     dofill(t,list_out_firing_rate,'blue',1,smval);
                    %                     xlabel('Time from Fixation (ms)')
                    %                     ylabel('Firing Rate')
                    %                     legend('In','Out','Location','Northwest')
                    %                     title('List Only')
                    %                     yl = ylim;
                    %                     if yl(1) < 0
                    %                         yl(1) = 0;
                    %                         ylim([0 yl(2)])
                    %                     end
                    %                     gaps = findgaps(find(list_sig_times{1,unit}));
                    %                     if ~isempty(gaps)
                    %                         for g = 1:size(gaps,1)
                    %                             gp = gaps(g,:);
                    %                             gp(gp == 0) = [];
                    %                             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    %                                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                    %                             uistack(h,'down')
                    %                             set(h,'facealpha',.25,'EdgeColor','None')
                    %                         end
                    %                     end
                    %                     hold off
                    %
                    %                     subplot(2,3,2)
                    %                     dofill(t,sequence_firing_rate,'red',1,smval);
                    %                     hold on
                    %                     dofill(t,out_sequence_firing_rate,'blue',1,smval);
                    %                     xlabel('Time from Fixation (ms)')
                    %                     ylabel('Firing Rate')
                    %                     legend('In','Out','Location','Northwest')
                    %                     title('Sequence Only')
                    %                     yl = ylim;
                    %                     if yl(1) < 0
                    %                         yl(1) = 0;
                    %                         ylim([0 yl(2)])
                    %                     end
                    %                     gaps = findgaps(find(nss));
                    %                     if ~isempty(gaps)
                    %                         for g = 1:size(gaps,1)
                    %                             gp = gaps(g,:);
                    %                             gp(gp == 0) = [];
                    %                             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    %                                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                    %                             uistack(h,'down')
                    %                             set(h,'facealpha',.25,'EdgeColor','None')
                    %                         end
                    %                     end
                    %                     hold off
                    %
                    %
                    %                     subplot(2,3,3)
                    %                     dofill(t,sequence_firing_rate,'green',1,smval);
                    %                     hold on
                    %                     dofill(t,list_in_firing_rate,'black',1,smval);
                    %                     xlabel('Time from Fixation (ms)')
                    %                     ylabel('Firing Rate')
                    %                     legend('Seq','List','Location','Northwest')
                    %                     title('List vs Seq')
                    %                     yl = ylim;
                    %                     if yl(1) < 0
                    %                         ylim([0 yl(2)])
                    %                     end
                    %                     gaps = findgaps(find(task_sig_times));
                    %                     if ~isempty(gaps)
                    %                         for g = 1:size(gaps,1)
                    %                             gp = gaps(g,:);
                    %                             gp(gp == 0) = [];
                    %                             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    %                                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                    %                             uistack(h,'down')
                    %                             set(h,'facealpha',.25,'EdgeColor','None')
                    %                         end
                    %                     end
                    %                     hold off
                    %                     title('List In vs Sequence In')
                    %
                    %                     subplot(2,3,4)
                    %                     [nl,nnl] = hist(max_list,20);
                    %                     [ns,nns] = hist(max_seq,nnl);
                    %                     hold on
                    %                     plot(nnl,filtfilt(1/3*ones(1,3),1,nl/sum(nl)),'k')
                    %                     plot(nns,filtfilt(1/3*ones(1,3),1,ns/sum(ns)),'g')
                    %                     hold off
                    %                     xlabel('Firing Rate')
                    %                     ylabel('Trial Count')
                    %                     title(['Peak FR diff: ' num2str(100*sum(abs(max_diff) > abs(shuff_max_diff))/numshuffs,3) '%']);
                    %
                    %                     subplot(2,3,5)
                    %                     [nl,nnl] = hist(window_list,20);
                    %                     [ns,nns] = hist(window_seq,nnl);
                    %                     hold on
                    %                     plot(nnl,filtfilt(1/3*ones(1,3),1,nl/sum(nl)),'k')
                    %                     plot(nns,filtfilt(1/3*ones(1,3),1,ns/sum(ns)),'g')
                    %                     hold off
                    %                     xlabel('Firing Rate')
                    %                     ylabel('Trial Count')
                    %                     title(['Window FR diff: ' num2str(100*sum(abs(window_diff) > abs(shuff_window_diff))/numshuffs,3) '%']);
                    %
                    %                     if sum(abs(max_diff) > abs(shuff_max_diff))/numshuffs > 0.95
                    %                         fr_task_context_count1 = [fr_task_context_count1 1];
                    %                     else
                    %                         fr_task_context_count1 = [fr_task_context_count1 0];
                    %                     end
                    %
                    %                     if sum(abs(window_diff) > abs(shuff_window_diff))/numshuffs > 0.95
                    %                         fr_task_context_count2 = [fr_task_context_count2 1];
                    %                         if window_diff > 0; %prefers list task
                    %                             fr_context_faster_slower = [fr_context_faster_slower 1];
                    %                         else %prefers sequence task
                    %                             fr_context_faster_slower = [fr_context_faster_slower -1];
                    %                         end
                    %                     else
                    %                         fr_task_context_count2 = [fr_task_context_count2 0];
                    %                         fr_context_faster_slower = [fr_context_faster_slower 0];
                    %                     end
                    %
                    %                     substring = [task_file(1:end-11) ' ' unit_stats{1,unit}];
                    %                     if sig_p_seq(end) == 1
                    %                         substring = [substring ' spatially consistent across tasks'];
                    %                     elseif  sig_p_seq(end) == -1
                    %                            substring = [substring ' spatially inconsistent across tasks'];
                    %                     else
                    %                         substring = [substring ' no spatial tuning in sequence task'];
                    %                     end
                    %                     subtitle(substring);
                    
                    %%
                    %save_and_close_fig('C:\Users\seth.koenig\Desktop\Significant Units\',[task_file(1:end-11) '_' unit_stats{1,unit}]);
                else
                    task_context_sig_count = [task_context_sig_count NaN];
                    fr_task_context_count1 = [fr_task_context_count1 NaN];
                    fr_task_context_count2 = [fr_task_context_count2 NaN];
                    fr_context_faster_slower = [fr_context_faster_slower NaN];
                end
            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                
                %---Misc Parameters---%
                monkey_all_unit_count(2,monk) = monkey_all_unit_count(2,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                
                %---Spatial Correlation Values for Non-Place cells---%
                all_non_place_cell_spatial_corrs = [all_non_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_non_place_cell_skagg_percents = [all_non_place_cell_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                non_place_skagg_names = [non_place_skagg_names  {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                
                subregion = session_data{sess}.subregion;
                vals = textscan(subregion,'%s','Delimiter',',');
                if length(vals{1}) == 1
                    non_place_cell_subregion = [non_place_cell_subregion vals{1}];
                else
                    chan = str2double(unit_stats{1,unit}(6));
                    non_place_cell_subregion = [non_place_cell_subregion vals{1}(chan)];
                end
            end
        end
    end
end
%% Copy Place Cell Figures to Summary Direction
clc

%---First Display Summary Results---%
disp(['Found ' num2str(length(all_list_peak_times)) ' place cells']);%number of place cells
disp(['Found ' num2str(sum(sig_p_list == 0)) ' place cells do not show significant increase in firing rate in field...removing'])
nans = find(sig_p_list == 0);
for i = 1:length(nans)
    disp(['Removing ' all_place_cell_unit_names{nans(i)}])
end
disp(['Found ' num2str(sum(isnan(all_list_peak_times))) ' place cells do not show a peak in time...removing'])
nans = find(isnan(all_list_peak_times));
for i = 1:length(nans)
    disp(['Removing ' all_place_cell_unit_names{nans(i)}])
end
%% Remove Un-Reliable Neurons and Neurons wthout a definitve peak in firing rate

%remove neurons without a definitive peak
all_context_gain(isnan(all_list_peak_times)) = [];
all_context_gain2(isnan(all_list_peak_times)) = [];
all_peak_seq(isnan(all_list_peak_times)) = [];
all_seq_peak_times(isnan(all_list_peak_times)) = [];
sig_p_seq(isnan(all_list_peak_times)) = [];
sig_p_list(isnan(all_list_peak_times)) = [];
all_peak_list(isnan(all_list_peak_times)) = [];
task_context_sig_count(isnan(all_list_peak_times)) = [];
%fr_task_context_count1(isnan(all_list_peak_times)) = [];
%fr_task_context_count2(isnan(all_list_peak_times)) = [];
%fr_context_faster_slower(isnan(all_list_peak_times)) = [];
all_whole_session_firing_rate(isnan(all_list_peak_times)) = [];
putative_excite_inhibit(isnan(all_list_peak_times)) = [];
all_seq_firing_curves(isnan(all_list_peak_times),:) = [];
all_in_rates(isnan(all_list_peak_times),:) = [];
all_list_peak_times(isnan(all_list_peak_times)) = [];
%%
tme = -twin1:twin2-1;
figure
subplot(2,2,1)
plot(tme,nanmean(all_in_rates))
hold on
plot(tme,nanmean(all_seq_firing_curves(sig_p_seq == ~0,:)))
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
plot([-44 -44],[yl(1) yl(2)],'k--')
plot([-twin1 twin2],[0 0],'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title(['Average View Cell Fixation Aligned Firing Rate'])
legend('Images Out2in','Seq Out2In')
xlim([-twin1 twin2])
box off
axis square

asfc = all_seq_firing_curves(sig_p_seq == ~0,:);
asfc(isnan(asfc(:,1)),:) = [];
subplot(2,2,2)
vals = asfc;
[~,mxi] = max(asfc');
[~,place_order] = sort(mxi); %sort order by peak firing time
imagesc([-twin1:twin2-1],[1:size(asfc,1)],asfc(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(asfc,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

%%

listFWHM = NaN(1,size(all_in_rates,1));
for n = 1:size(all_in_rates,1)
    cv = all_in_rates(n,:);
    peak = all_list_peak_times(n);
    tims = find(cv < 0.5*cv(peak));
    diffs = tims-peak;
    neg_diffs = diffs(diffs < 0);
    pos_diffs = diffs(diffs > 0);
    neg = min(abs(neg_diffs));
    pos = min(pos_diffs);
    
    if isempty(pos)
        pos = twin2;
        listFWHM(n) = NaN;
    else
        listFWHM(n) = neg+pos;
    end
    
end
%%
asfc = all_seq_firing_curves(sig_p_seq == ~0,:);
aspt = all_seq_peak_times(sig_p_seq == ~0);
aspt(isnan(asfc(:,1))) = [];
asfc(isnan(asfc(:,1)),:) = [];

seqFWHM = NaN(1,size(asfc,1));
for n = 1:size(asfc,1)
    cv = asfc(n,:);
    if ~isnan(aspt(n))
        peak = aspt(n);
    else
        [~,peak] = max(cv);
    end
    tims = find(cv < 0.5*cv(peak));
    diffs = tims-peak;
    neg_diffs = diffs(diffs < 0);
    pos_diffs = diffs(diffs > 0);
    neg = min(abs(neg_diffs));
    pos = min(pos_diffs);
    
    if isempty(pos)
        pos = twin2;
        seqFWHM(n) = NaN;
    elseif isempty(neg);
        seqFWHM(n) = NaN;
    else
        seqFWHM(n) = neg+pos;
    end
end
%%
subplot(2,2,3)
histogram(listFWHM,25,'Normalization','probability')
hold on
histogram(seqFWHM,25,'Normalization','probability')
hold off
xlabel('FWHM (ms)')
ylabel('Proportion')
legend('List','Seq')