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
num_shuffs_crossvalid = 1000;
task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800;
imageY = 600;


%---Fixation/Saccade Duration parameters---%
min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

only_process_active_units = false; %units with criterion below
exclude_nonplace_eye_cells = false; %exclude saccade modulated units that are not place cells as controls
t_start = 50; %time before fixation to include in information score calculation, this period is not warped
min_saccades_with_spikes = 0.05;% 5% to remove neurons that have basically no activity since they pass too :(...
%essentially 0.4 Hz threshold if window is 50 ms in wide
min_num_fixations = 100;
min_num_fixations_w_spikes = 100;

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
        
        %load ege movement modulation analysis data
        load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
        
        %should also load direction modulation
        load([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'])
        
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
                if only_process_active_units
                    %                 if sum(spike_counts_out2in > 0) < min_saccades_with_spikes*size(firing_rate,1) || ...
                    %                         min_num_fixations > size(firing_rate,1)
                    %
                    %                         skipped_units = skipped_units+1;
                    %                         continue
                    %                     end
                    %                 end
                    if sum(spike_counts_out2in > 0) < min_num_fixations_w_spikes
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
                [~,all_out2in{place_cell_ind}] = nandens(firing_rate,2*smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
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
                
                
                if exclude_nonplace_eye_cells
                    %remove other types of eye movement modulated units
                    if mrls.all_saccades_shuffled_prctile(unit) > 95 || amplitude_correlations_percentile(unit) > 97.5
                        skipped_units = skipped_units+1;
                        continue
                    elseif (temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) > 95) ... %significant stability
                            && (temporal_info.fixation.shuffled_rate_prctile(unit) > 95)
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
                [~,all_out2in_nonplace{nonplace_cell_ind}] = nandens(firing_rate,2*smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                nonplace_cell_ind = nonplace_cell_ind+1;
            end
        end
    end
end


%do some cleanup
clc
all_out2in(place_cell_ind:end) = [];
all_out2in_nonplace(nonplace_cell_ind:end) = [];

%%
%% Run Shuffling Methods

crossValidatedStatsPlace = crossValidatePlacePopulationTemporalResponse(all_out2in,num_shuffs_crossvalid,'Place',twin1,twin2);
crossValidatedStatsNonPlace = crossValidatePlacePopulationTemporalResponse(all_out2in_nonplace,num_shuffs_crossvalid,'NonPlace',twin1,twin2);

