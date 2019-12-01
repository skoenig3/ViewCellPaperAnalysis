%Code Written to Make Time Cell Plots Limited by Fixation Duration Rather
%than blindly plotting from fixation start to twin (usually 400 ms).
%Code Is essentially the same as Population_Place_Cell_Fixation_Anlaysis
%otherwise.
%Written by Seth Koenig 5/15/19

clar %clear,clc

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s

%---Spatial Correlation Values for all cells---%
all_place_cell_spatial_corrs = []; %place cell spatial correlations
all_place_cell_spatial_skagg_percents= []; %place cell skagg information score percentiles
all_non_place_cell_spatial_corrs = [];%non-place cell spatial correlations
all_non_place_cell_skagg_percents = [];%non-place cell skagg information score percentiles
non_place_skagg_names = [];

%---Place Cell Firing Rate Curves for List---%
%first, median, mean
all_out2in_spike_times = [];
all_out2out_spike_times = [];
all_in2out_spike_times = [];
all_in2in_spike_times = [];
all_in2in_first_spike_times = [];
all_in2in_second_spike_times = [];

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
        
        %load spatial analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        %Save as new variableso can reload later...kind of dumb but that was how it was written
        absolute_fixationstats = fixationstats;
        absolute_cfg = cfg;
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{sess},task);
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        %remove units with too few trials
        %these are the absolute minimum data required to do data analysis
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if isnan(stats_across_tasks(1,unit))%no peak detected
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                
                %---Misc Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s for monkey
                
                %---Spatial Correlation Values for Place cells---%
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_place_cell_spatial_skagg_percents = [all_place_cell_spatial_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                
                
                place_field_matrix = all_place_field_matrix{unit};
                %note ignores incomplete coverage
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Calculate Firing Rate Locked to Fixations---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fix_in_out = NaN(1,1000); %firing rate for fixations inside or outside  of place field
                %1) first fixation in: out-> in
                %2) fixation in field but not first: in -> in
                %3) first fixation out of field: in -> out
                %4) fixation out of field but not first: out-> out
                
                next_fix_in_out = NaN(1,1000);
                all_fix_dur = NaN(1,1000);
                all_previous_fix_dur = NaN(1,1000);%include prior sacdur
                fix_ind = 1; %fixation # so can track in variables above
                
                fixationstats = absolute_fixationstats; %reload because written over below
                cfg = absolute_cfg; %reload because written over below
                num_trials = length(cfg.trl);%number of trials
                
                
                out2in_spike_times = NaN(1000,3);
                out2out_spike_times = NaN(1000,3);
                in2out_spike_times = NaN(1000,3);
                in2in_spike_times = NaN(1000,3);
                in2in_first_spike_times = NaN(1000,3);
                in2in_second_spike_times = NaN(1000,3);

                for t = 1:num_trials
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                        if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                            imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                            
                            %---image info---%
                            img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                            
                            % if monkey isn't paying attention and looked away image presentation
                            % is now longer than imgdur (because of cumulative looking time)
                            % so data isn't probably worth much plus have to cut off somewhere
                            if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                                imgoff = imgon+1.5*imgdur-1;
                            end
                            imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                            
                            fixationtimes = fixationstats{t}.fixationtimes; %fixation start and end times
                            saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                            fixations = round(fixationstats{t}.fixations); %mean fixation location
                            xy = fixationstats{t}.XY; %xy eye trace
                            
                            %find fiations and saccades that did not occur during the image period;
                            %should also take care of the 1st fixation on the crosshair
                            
                            %fixation started before image turned on
                            invalid= find(fixationtimes(1,:) < imgon);
                            fixationtimes(:,invalid) = [];
                            fixations(:,invalid) = [];
                            
                            %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                            invalid= find(fixationtimes(1,:) > imgoff-twin2);
                            fixationtimes(:,invalid) = [];
                            fixations(:,invalid) = [];
                            
                            %saccade started before image turned on
                            invalid= find(saccadetimes(1,:) < imgon);
                            saccadetimes(:,invalid) = [];
                            
                            %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                            invalid= find(saccadetimes(1,:) > imgoff-twin2);
                            saccadetimes(:,invalid) = [];
                            
                            %remove fixations that are too short to really look at for analysis
                            fixdur = fixationtimes(2,:)-fixationtimes(1,:)+1;
                            fixationtimes(:,fixdur < min_fix_dur) = [];
                            fixations(:,fixdur < min_fix_dur) = [];
                            
                            img_index = find(img_cnd == cfg.trl(t).cnd);
                            %which image monkey viewed if nan or empty then skip cuz
                            %bad trial
                            if isempty(img_index) || any(isnan(which_img(img_index)))
                                continue
                            end
                            
                            spikes = find(data(unit).values{t}); %spike trains for this trial
                            for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                                fixdur = fixationtimes(2,f)-fixationtimes(1,f);
                                prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                                if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                                    continue;
                                end
                                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                                if sacamp < min_sac_amp %prior saccade is too small so ignore
                                    continue
                                end
                                                                
                                
                                prior_fix_in_out = NaN;
                                %determine if prior fixation was in or out of place field
                                fixx = fixations(1,f-1);
                                fixx(fixx < 1) = 1;
                                fixx(fixx > imageX) = imageX;
                                fixy = imageY-fixations(2,f-1);
                                fixy(fixy < 1) = 1;
                                fixy(fixy > imageY) = imageY;
                                
                                if place_field_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                                    prior_fix_in_out = 1;
                                else
                                    prior_fix_in_out = 0;
                                end
                                last_fixx = fixx;
                                last_fixy = fixy;
                                
                                %determine if next fixation was in or out of place field
                                if f~= size(fixations,2)
                                    fixx = fixations(1,f+1);
                                    fixx(fixx < 1) = 1;
                                    fixx(fixx > imageX) = imageX;
                                    fixy = imageY-fixations(2,f+1);
                                    fixy(fixy < 1) = 1;
                                    fixy(fixy > imageY) = imageY;
                                    
                                    if place_field_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                                        next_fix_in_out(fix_ind) = 1;
                                    else
                                        next_fix_in_out(fix_ind) = 0;
                                    end
                                    last_fixx = fixx;
                                    last_fixy = fixy;
                                end
                                
                                %determine if this fixation was in our out of place field
                                fixx = fixations(1,f);
                                fixx(fixx < 1) = 1;
                                fixx(fixx > imageX) = imageX;
                                fixy = imageY-fixations(2,f);
                                fixy(fixy < 1) = 1;
                                fixy(fixy > imageY) = imageY;
                                first_in2in = false;
                                if place_field_matrix(fixy,fixx) == 1 %then inside
                                    if prior_fix_in_out == 1%prior fixation was inside so in->in
                                        fix_in_out(fix_ind) = 2;
                                        if f > 2 && fix_ind > 1
                                            first_in2in = true;
                                        end
                                    else %out->in
                                        fix_in_out(fix_ind) = 1;
                                    end
                                elseif place_field_matrix(fixy,fixx) == 0 %then inside, NaNs are for locations not analyzed
                                    if prior_fix_in_out == 1%prior fixation was inside so in->out
                                        fix_in_out(fix_ind) = 3;
                                    else %out->out
                                        fix_in_out(fix_ind) = 4;
                                    end
                                else %not a valid fixation location too sparse of coverage to analyze
                                    continue
                                end
                                
                                
                                fixt = fixationtimes(1,f);%start of fixation
                                fix_spikes = spikes(spikes > fixt & spikes <= fixt+fixdur)-fixt;
                                
                                if fix_in_out(fix_ind) == 1 %out2in
                                    if ~isempty(fix_spikes)
                                        out2in_spike_times(fix_ind,1) = fix_spikes(1);
                                        out2in_spike_times(fix_ind,2) = median(fix_spikes);
                                        out2in_spike_times(fix_ind,3) = mean(fix_spikes);
                                    end
                                elseif fix_in_out(fix_ind) == 2%in2in
                                    if ~isempty(fix_spikes)
                                        in2in_spike_times(fix_ind,1) = fix_spikes(1);
                                        in2in_spike_times(fix_ind,2) = median(fix_spikes);
                                        in2in_spike_times(fix_ind,3) = mean(fix_spikes);
                                        
                                        if first_in2in
                                            if fix_in_out(fix_ind-1) == 1 %prior was out2in
                                                in2in_first_spike_times(fix_ind,1) = fix_spikes(1);
                                                in2in_first_spike_times(fix_ind,2) = median(fix_spikes);
                                                in2in_first_spike_times(fix_ind,3) = mean(fix_spikes);
                                            else
                                                in2in_second_spike_times(fix_ind,1) = fix_spikes(1);
                                                in2in_second_spike_times(fix_ind,2) = median(fix_spikes);
                                                in2in_second_spike_times(fix_ind,3) = mean(fix_spikes);
                                            end
                                        end
                                    end
                                elseif fix_in_out(fix_ind) == 3%out2in
                                    if ~isempty(fix_spikes)
                                        in2out_spike_times(fix_ind,1) = fix_spikes(1);
                                        in2out_spike_times(fix_ind,2) = median(fix_spikes);
                                        in2out_spike_times(fix_ind,3) = mean(fix_spikes);
                                    end
                                elseif fix_in_out(fix_ind) == 4 %out2out
                                    if ~isempty(fix_spikes)
                                        out2out_spike_times(fix_ind,1) = fix_spikes(1);
                                        out2out_spike_times(fix_ind,2) = median(fix_spikes);
                                        out2out_spike_times(fix_ind,3) = mean(fix_spikes);
                                    end
                                end
                                
                                fix_ind = fix_ind+1;
                            end
                        end
                    end
                end
                
                all_out2in_spike_times = [all_out2in_spike_times; nanmean(out2in_spike_times)];
                all_out2out_spike_times = [all_out2out_spike_times; nanmean(out2out_spike_times)];
                all_in2out_spike_times = [all_in2out_spike_times; nanmean(in2out_spike_times)];
                all_in2in_spike_times = [all_in2in_spike_times; nanmean(in2in_spike_times)];
                
                all_in2in_first_spike_times = [all_in2in_first_spike_times; nanmean(in2in_first_spike_times)];
                all_in2in_second_spike_times = [all_in2in_second_spike_times; nanmean(in2in_second_spike_times)];
            end
        end
    end
end
%%

avg_data = [nanmean(all_out2in_spike_times); nanmean(all_in2in_first_spike_times); ...
    nanmean(all_in2in_second_spike_times); nanmean(all_in2in_spike_times); ...
    nanmean(all_in2out_spike_times); nanmean(all_out2out_spike_times)];

figure
bar(avg_data')
xlabel('Method')
set(gca,'xticklabels',{'First','Median','Mean'})
ylabel('Spike Latency (ms)')
legend({'Out2in','In2In 1st','In2In 2nd+','In2In All','In2Out','Out2Out'})
box off
title('View Cell Population Spike Latencys')

%%
num_cells = size(all_out2in_spike_times,1);
outin_grouping = [ones(num_cells,1); 2*ones(num_cells,1); 3*ones(num_cells,1); ...
    4*ones(num_cells,1); 5*ones(num_cells,1); 6*ones(num_cells,1)];

%first spike
[P1,~,STATS1] = anova1([all_out2in_spike_times(:,1); all_in2in_first_spike_times(:,1); ...
    all_in2in_second_spike_times(:,1); all_in2in_spike_times(:,1); ...
    all_in2out_spike_times(:,1); all_out2out_spike_times(:,1)],outin_grouping);
COMPARISON1 = multcompare(STATS1);
%%
%median spike
[P2,~,STATS2] = anova1([all_out2in_spike_times(:,2); all_in2in_first_spike_times(:,2); ...
    all_in2in_second_spike_times(:,2); all_in2in_spike_times(:,2); ...
    all_in2out_spike_times(:,2); all_out2out_spike_times(:,2)],outin_grouping);
COMPARISON2 = multcompare(STATS2);
%%
%mean spike
[P3,~,STATS3] = anova1([all_out2in_spike_times(:,3); all_in2in_first_spike_times(:,3); ...
    all_in2in_second_spike_times(:,3); all_in2in_spike_times(:,3); ...
    all_in2out_spike_times(:,3); all_out2out_spike_times(:,3)],outin_grouping);
COMPARISON3 = multcompare(STATS3);
    