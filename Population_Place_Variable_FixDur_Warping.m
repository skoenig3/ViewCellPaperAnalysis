%Code cleaned up and re-checked, 2/17/20 after many modification suggestions
%noted a few minor indexing errors in LFP calculations
%Code modified from Draft Version in ListSQ on 6/3/19
%Code determines if spatiotemporal responses are better aligned to
%fixation/saccade onset in a linear manner or are scaled by
%fixation/saccade duration. Code also anlayzes LFPs in this way too.

clar %clear,clc

%---Figure Options---%
plot_individual_cells = true;
set(0,'DefaultFigureVisible','ON');
figure_dir2 = 'C:\Users\seth.koenig\Desktop\New folder\';

%---Import Options---%
task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800;
imageY = 600;
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial

Fs = 1000; %Hz sampling frequency
fixwin = 5;%size of fixation window on each crosshair
smval = 30;%2*std of gaussian kernel so 15 ms standard deviation
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

%---Fixation/Saccade Duration parameters---%
cutoff_short_fixations = true; %cutt off spikes after fixation has ended 

min_saccades_with_spikes = 0.05;% 5% to remove neurons that have basically no activity since they pass too :(...
%essentially 0.4 Hz threshold if window is 50 ms in wide
min_num_fixations = 100;

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

t_start = 50; %time before fixation to include in information score calculation, this period is not warped
max_fix_dur = 600;%max fixation duration to warp, shouldn't be to many but could include 2nd minisaccades
fix_dur_to_warp_to = 180;%1 ms/degree of phase @ 180 ms, close to median fix dur of 188


%---Storage Variables---%
%resampling expected to lead to no change in firing rate, while scaling
%make cause changes in firing rate due to longtail fixation duration distribution
%Unwarped and Warped x out2in and all fixations
%row 1, ou2in raw
%row 2, out2in warped
%row 3, all raw
%row 4, all warped
fixation_skagg_info_raw_resampled = NaN(4,110);%fixation linear-spaced up & down sampling
fixation_LFP_skagg_info_raw_resampled = NaN(2,350);%LFP linear-spaced up & down sampling
fixation_skagg_info_raw_scaled = NaN(4,110); %fixation time stamp scaling


%determine if firing rate or amplitude has changed with scaling organized as above
%for before/after warping
fixation_average_firing_rate_resampled = NaN(4,110);%firing rates for resampling
fixation_LFP_amplitude_resampled = NaN(2,350);%LFP amplitudes for resampling
fixation_average_firing_rate_scaled = NaN(4,110);%firing rates for scaling



%----Main Script---%
%timing variables that are the same for all
t11 = -twin1:fix_dur_to_warp_to-1;
total_time = 1:fix_dur_to_warp_to+t_start+1;
p_x = total_time/sum(total_time);

cell_ind = 1;
LFP_ind = 1;
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
        

        %load Place Cell Fixation Analysis data
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
        
        for unit =1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            if isnan(stats_across_tasks(1,unit))%no peak detected
                continue
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Process Spike Data---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                place_field_matrix = all_place_field_matrix{unit};
                %note ignores incomplete coverage
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Calculate Firing Rate Locked to Fixations---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
                %1) first fixation in: out-> in
                %2) fixation in field but not first: in -> in
                %3) first fixation out of field: in -> out
                %4) fixation out of field but not first: out-> out
                
                %out2in fixations only
                out2in_fix_dur = NaN(1,1000);
                fix_locked_firing_raw = NaN(1000,(twin1+fix_dur_to_warp_to)); %spike trains locked to fixations, unwarped
                fix_locked_firing_resampled = NaN(1000,(twin1+fix_dur_to_warp_to)); %resampled spike trains locked to fixations
                fix_locked_firing_scaled = NaN(1000,(twin1+fix_dur_to_warp_to)); %scaled spike trains locked to fixations
                
                %all fixations
                all_fix_dur = NaN(1,5000);
                all_fix_locked_firing_raw = NaN(5000,(twin1+fix_dur_to_warp_to)); %all spike trains locked to fixations
                all_fix_locked_firing_resampled = NaN(5000,(twin1+fix_dur_to_warp_to));%resampled spike trains locked to fixations
                all_fix_locked_firing_scaled = NaN(5000,(twin1+fix_dur_to_warp_to));  %scaled spike trains locked to fixations
                
                
                fix_ind = 1; %fixation # so can track in variables above
                all_fix_ind = 1; %fixation # so can track in variables above
                
                fixationstats = absolute_fixationstats; %reload because written over below
                cfg = absolute_cfg; %reload because written over below
                num_trials = length(cfg.trl);%number of trials
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
                            fixdur(fixdur < min_fix_dur) = [];
                            
                            fixationtimes(:,fixdur > max_fix_dur) = [];
                            fixations(:,fixdur > max_fix_dur) = [];
                            fixdur(fixdur > max_fix_dur) = [];
                            
                            img_index = find(img_cnd == cfg.trl(t).cnd);
                            %which image monkey viewed if nan or empty then skip cuz
                            %bad trial
                            if isempty(img_index) || any(isnan(which_img(img_index)))
                                continue
                            end
                            
                            spikes = find(data(unit).values{t}); %spike trains for this trial
                            for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                                fixt = fixationtimes(1,f);%start of fixation
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
                                
                                %determine if this fixation was in our out of place field
                                fixx = fixations(1,f);
                                fixx(fixx < 1) = 1;
                                fixx(fixx > imageX) = imageX;
                                fixy = imageY-fixations(2,f);
                                fixy(fixy < 1) = 1;
                                fixy(fixy > imageY) = imageY;
                                if place_field_matrix(fixy,fixx) == 1 %then inside
                                    if prior_fix_in_out == 1%prior fixation was inside so in->in
                                        fix_in_out(fix_ind) = 2;
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
                                
                                %---For Out2In Fixations Only---%
                                if fix_in_out(fix_ind) == 1 %out2in
                                    %clear temp varaibles just in case...should've renamed these
                                    temp = [];
                                    temp2 = [];
                                    temp3 = [];
                                    
                                    %get firing rate locked to fixation
                                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+fix_dur_to_warp_to)-fixt+twin1;
                                    temp = zeros(1,twin1+fix_dur_to_warp_to);
                                    temp(fix_spikes) = 1;
                                    if fixdur < fix_dur_to_warp_to && cutoff_short_fixations
                                        remaining = fix_dur_to_warp_to-fixdur;
                                        temp(end-remaining+1:end) = NaN;
                                    end
                                    fix_locked_firing_raw(fix_ind,:) = temp;
                                    
                                    if fixdur == fix_dur_to_warp_to
                                        fix_locked_firing_resampled(fix_ind,:) = temp;
                                        fix_locked_firing_scaled(fix_ind,:) = temp;
                                    else
                                        fix_locked_firing_resampled(fix_ind,1:twin1) = temp(1:twin1);
                                        fix_locked_firing_scaled(fix_ind,1:twin1) = temp(1:twin1);
                                        
                                        %resampling
                                        fix_spikes2 = spikes(spikes > fixt & spikes <= fixt + fixdur)-fixt;
                                        temp2 = zeros(1,fixdur);
                                        temp2(fix_spikes2) = 1;
                                        resampled = round(linspace(1,fixdur,fix_dur_to_warp_to));
                                        temp2 = temp2(resampled);
                                        fix_locked_firing_resampled(fix_ind,twin1+1:end) = temp2;
                                        
                                        %scaling
                                        fix_spikes2_scaled = round(fix_spikes2*fix_dur_to_warp_to/fixdur);
                                        fix_spikes2_scaled(fix_spikes2_scaled == 0) = [];
                                        fix_spikes2_scaled(fix_spikes2_scaled > fix_dur_to_warp_to) = [];
                                        temp3 = zeros(1,fix_dur_to_warp_to);
                                        temp3(fix_spikes2_scaled) = 1;
                                        fix_locked_firing_scaled(fix_ind,twin1+1:end) = temp3;
                                    end
                                    out2in_fix_dur(fix_ind) = fixdur;
                                    fix_ind = fix_ind+1;
                                end
                                
                                %---For All Fixations---%
                                %clear temp varaibles just in case...should've renamed these
                                temp = [];
                                temp2 = [];
                                temp3 = [];
                                fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+fix_dur_to_warp_to)-fixt+twin1;
                                temp = zeros(1,twin1+fix_dur_to_warp_to);
                                temp(fix_spikes) = 1;
                                if fixdur < fix_dur_to_warp_to && cutoff_short_fixations
                                    remaining = fix_dur_to_warp_to-fixdur;
                                    temp(end-remaining+1:end) = NaN;
                                end
                                all_fix_locked_firing_raw(all_fix_ind,:) = temp;
                                
                                if fixdur == fix_dur_to_warp_to
                                    all_fix_locked_firing_resampled(all_fix_ind,:) = temp;
                                    all_fix_locked_firing_scaled(all_fix_ind,:) = temp;
                                else
                                    all_fix_locked_firing_resampled(all_fix_ind,1:twin1) = temp(1:twin1);
                                    all_fix_locked_firing_scaled(all_fix_ind,1:twin1) = temp(1:twin1);
                                    
                                    %resampling
                                    fix_spikes2 = spikes(spikes > fixt & spikes <= fixt + fixdur)-fixt;
                                    temp2 = zeros(1,fixdur);
                                    temp2(fix_spikes2) = 1;
                                    resampled = round(linspace(1,fixdur,fix_dur_to_warp_to));
                                    temp2 = temp2(resampled);
                                    all_fix_locked_firing_resampled(all_fix_ind,twin1+1:end) = temp2;
                                    
                                    %scaling
                                    all_fix_spikes2_scaled = round(fix_spikes2*fix_dur_to_warp_to/fixdur);
                                    all_fix_spikes2_scaled(all_fix_spikes2_scaled == 0) = [];
                                    all_fix_spikes2_scaled(all_fix_spikes2_scaled > fix_dur_to_warp_to) = [];
                                    temp3 = zeros(1,fix_dur_to_warp_to);
                                    temp3(all_fix_spikes2_scaled) = 1;
                                    all_fix_locked_firing_scaled(all_fix_ind,twin1+1:end) = temp3;
                                    
                                    if abs(sum(temp3)-length(fix_spikes2)) > 2 %2 for rounding cuz could loose 1st and last spikes
                                        disp('now')
                                    end
                                end
                                all_fix_dur(all_fix_ind) = fixdur;
                                all_fix_ind = all_fix_ind+1; %fixation # so can track in variables above
                            end
                        end
                    end
                end
                
                %---Remove excess NaNs---%
                out2in_fix_dur = laundry(out2in_fix_dur);
                fix_locked_firing_raw = laundry(fix_locked_firing_raw);
                fix_locked_firing_resampled = laundry(fix_locked_firing_resampled);
                fix_locked_firing_scaled = laundry(fix_locked_firing_scaled);
                                
                all_fix_locked_firing_raw = laundry(all_fix_locked_firing_raw);
                all_fix_locked_firing_resampled = laundry(all_fix_locked_firing_resampled);
                all_fix_locked_firing_scaled = laundry(all_fix_locked_firing_scaled);
                all_fix_dur = laundry(all_fix_dur);
                
                
                %--Determine if should NOT process Data do to sparse sampling---%
                spike_counts_out2in = nansum(fix_locked_firing_raw(:,twin1-t_start+1:end));
                if sum(spike_counts_out2in > 0) < min_saccades_with_spikes*size(fix_locked_firing_raw,1) || ...
                        min_num_fixations > size(fix_locked_firing_raw,1)
                    process_out2in = false;
                else
                    process_out2in = true;
                end
                
                spike_counts_all = nansum(all_fix_locked_firing_raw(:,twin1-t_start+1:end));
                if sum(spike_counts_all > 0) < min_saccades_with_spikes*size(all_fix_locked_firing_raw,1) || ...
                        min_num_fixations > size(all_fix_locked_firing_raw,1)
                    process_all = false;
                else
                    process_all = true;
                end
                  
                %---Calculate Temporal Information for Raw & Warped Data---%
                if  process_out2in
                    if cutoff_short_fixations
                        yraw = nandens3(fix_locked_firing_raw(:,twin1-50:end),smval,Fs);
                        yresampled = nandens3(fix_locked_firing_resampled(:,twin1-50:end),smval,Fs);
                        yscaled = nandens3(fix_locked_firing_scaled(:,twin1-50:end),smval,Fs);
                    else
                        yraw = nandens(fix_locked_firing_raw(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                        yresampled = nandens(fix_locked_firing_resampled(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                        yscaled = nandens(fix_locked_firing_scaled(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                    end
                    
                    lambda_x = yraw;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    raw_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    lambda_x = yresampled;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    resampled_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    lambda_x = yscaled;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    scaled_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    fixation_skagg_info_raw_resampled(1,cell_ind) = raw_skaggs;
                    fixation_skagg_info_raw_resampled(2,cell_ind) = resampled_skaggs;
                    
                    fixation_skagg_info_raw_scaled(1,cell_ind) = raw_skaggs;
                    fixation_skagg_info_raw_scaled(2,cell_ind) = scaled_skaggs;
                    
                    fixation_average_firing_rate_resampled(1,cell_ind) = mean(yraw(twin1-t_start+1:end));
                    fixation_average_firing_rate_resampled(2,cell_ind) = mean(yresampled(twin1-t_start+1:end));
                    
                    fixation_average_firing_rate_scaled(1,cell_ind) = mean(yraw(twin1-t_start+1:end));
                    fixation_average_firing_rate_scaled(2,cell_ind) = mean(yscaled(twin1-t_start+1:end));
                end
                
                if  process_all
                    if cutoff_short_fixations
                        yraw = nandens3(all_fix_locked_firing_raw(:,twin1-50:end),smval,Fs);%'gauss',Fs,'nanflt');
                        yresampled = nandens3(all_fix_locked_firing_resampled(:,twin1-50:end),smval,Fs);%'gauss',Fs,'nanflt');
                        yscaled = nandens3(all_fix_locked_firing_scaled(:,twin1-50:end),smval,Fs);%'gauss',Fs,'nanflt');
                    else
                        yraw = nandens(all_fix_locked_firing_raw(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                        yresampled = nandens(all_fix_locked_firing_resampled(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                        yscaled = nandens(all_fix_locked_firing_scaled(:,twin1-50:end),smval,'gauss',Fs,'nanflt');
                    end
                    
                    lambda_x = yraw;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    raw_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    lambda_x = yresampled;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    resampled_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    lambda_x = yscaled;
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    scaled_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    
                    fixation_skagg_info_raw_resampled(3,cell_ind) = raw_skaggs;
                    fixation_skagg_info_raw_resampled(4,cell_ind) = resampled_skaggs;
                    
                    fixation_skagg_info_raw_scaled(3,cell_ind) = raw_skaggs;
                    fixation_skagg_info_raw_scaled(4,cell_ind) = scaled_skaggs;
                    
                    fixation_average_firing_rate_resampled(3,cell_ind) = mean(yraw(twin1-t_start+1:end));
                    fixation_average_firing_rate_resampled(4,cell_ind) = mean(yresampled(twin1-t_start+1:end));
                    
                    fixation_average_firing_rate_scaled(3,cell_ind) = mean(yraw(twin1-t_start+1:end));
                    fixation_average_firing_rate_scaled(4,cell_ind) = mean(yscaled(twin1-t_start+1:end));
                end

                %%%---Plot Indivisual Cells---%%%
                if plot_individual_cells
                    %%
                    figure
                    
                    [~,si] = sort(out2in_fix_dur);
                    fix_locked_firing_raw = fix_locked_firing_raw(si,:);
                    [trial,time] = find(fix_locked_firing_raw == 1);
                    fix_locked_firing_resampled = fix_locked_firing_resampled(si,:);
                    [trialw,timew] = find(fix_locked_firing_resampled == 1);
                    
                    fix_locked_firing_scaled = fix_locked_firing_scaled(si,:);
                    [trials,times] = find(fix_locked_firing_scaled == 1);
                    
                    subplot(2,2,1)
                    plot(time-twin1,trial,'.k')
                    hold on
                    plot(timew-twin1,trialw+max(trial),'.b')
                    plot(times-twin1,trials+max(trial)+max(trialw),'.r')
                    hold off
                    ylim([0 max(trialw) + max(trial) + max(trials)])
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Fixation #')
                    title('Rasters for Out2In Only')
                    box off
                    
                    subplot(2,2,2)
                    hold on
                    if cutoff_short_fixations
                        dofill2(t11,fix_locked_firing_raw,'black',1,smval);
                        dofill2(t11,fix_locked_firing_resampled,'blue',1,smval);
                        dofill2(t11,fix_locked_firing_scaled,'red',1,smval);
                    else
                        dofill(t11,fix_locked_firing_raw,'black',1,smval);
                        dofill(t11,fix_locked_firing_resampled,'blue',1,smval);
                        dofill(t11,fix_locked_firing_resampled,'sampled',1,smval);
                    end
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Fixation Start (ms)')
                    xlim([-twin1 twin1])
                    p_change_resample = 100*(fixation_skagg_info_raw_resampled(2,cell_ind)...
                        -fixation_skagg_info_raw_resampled(1,cell_ind))/fixation_skagg_info_raw_resampled(1,cell_ind);
                    p_change_scaled = 100*(fixation_skagg_info_raw_scaled(2,cell_ind)...
                        -fixation_skagg_info_raw_scaled(1,cell_ind))/fixation_skagg_info_raw_scaled(1,cell_ind);
                    title(['Out2In Only: \Delta Rescaled = ' num2str(p_change_resample,3) '% ' ...
                        '\Delta Scaled = ' num2str(p_change_scaled,3) '%'])
                    legend('Raw','Resampled','Scaled')
                    
                    [~,si] = sort(all_fix_dur);
                    all_fix_locked_firing_raw = all_fix_locked_firing_raw(si,:);
                    [trial,time] = find(all_fix_locked_firing_raw == 1);
                    all_fix_locked_firing_resampled = all_fix_locked_firing_resampled(si,:);
                    [trialw,timew] = find(all_fix_locked_firing_resampled == 1);
                    all_fix_locked_firing_scaled = all_fix_locked_firing_scaled(si,:);
                    [trials,times] = find(all_fix_locked_firing_scaled == 1);
                    
                    subplot(2,2,3)
                    plot(time-twin1,trial,'.k')
                    hold on
                    plot(timew-twin1,trialw+max(trial),'.b')
                    plot(times-twin1,trials+max(trial)+max(trialw),'.r')
                    hold off
                    ylim([0 max(trialw) + max(trial) + max(trials)])
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Fixation #')
                    title('Rasters for All Fixations, sorted by fixation duration')
                    box off
                    
                    subplot(2,2,4)
                    hold on
                    if cutoff_short_fixations
                        dofill2(t11,all_fix_locked_firing_raw,'black',1,smval);
                        dofill2(t11,all_fix_locked_firing_resampled,'blue',1,smval);
                        dofill2(t11,all_fix_locked_firing_scaled,'red',1,smval);
                    else
                        dofill(t11,all_fix_locked_firing_raw,'black',1,smval);
                        dofill(t11,all_fix_locked_firing_resampled,'blue',1,smval);
                        dofill(t11,all_fix_locked_firing_resampled,'red',1,smval);
                    end
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Fixation Start (ms)')
                    xlim([-twin1 twin1])
                    p_change_resample = 100*(fixation_skagg_info_raw_resampled(4,cell_ind)...
                        -fixation_skagg_info_raw_resampled(3,cell_ind))/fixation_skagg_info_raw_resampled(3,cell_ind);
                    p_change_scaled = 100*(fixation_skagg_info_raw_scaled(4,cell_ind)...
                        -fixation_skagg_info_raw_scaled(3,cell_ind))/fixation_skagg_info_raw_scaled(3,cell_ind);
                    title(['All Fixations: \Delta Rescaled = ' num2str(p_change_resample,3) '% ' ...
                        '\Delta Scaled = ' num2str(p_change_scaled,3) '%'])
                    
                    subtitle([task_file(1:end-11) '_' unit_stats{1,unit}])
                    %%
                    save_and_close_fig(figure_dir2,[task_file(1:end-11) '_' unit_stats{1,unit} '_view_cell_fixation_time_warping']);
                end
                cell_ind = cell_ind +1;
            else
                continue
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Process LFP Data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %seperate from spikes because importing "all" LPF channels
        %remove bad LFP channels
        LFPchannels = find_desired_channels(cfg,'LFP');
        bad_channels = [];
        for channel = 1:4
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        LFPchannels(bad_channels) = NaN;
        if all(isnan(LFPchannels))
            continue %no analyses to do
        end
        
        for chan = 1:length(LFPchannels)
            if ~isnan(LFPchannels(chan))
                %only do all fixations since no out2in
                fix_LFP = NaN(5000,(twin1+fix_dur_to_warp_to)); %spike trains locked to fixations
                fix_LFP_resampled = NaN(5000,(twin1+fix_dur_to_warp_to)); %spike trains locked to fixations
                all_LFP_fixdur = NaN(1,5000);
                fix_ind = 1; %fixation # so can track in variables above
                
                %calculate RMS on trial by trials since interferes ruins
                %amplitude measurements
                fix_LFP_RMS = NaN(1,5000);
                fix_LFP_RMS_resampled  = NaN(1,5000);
                
                fixationstats = absolute_fixationstats; %reload because written over below
                cfg = absolute_cfg; %reload because written over below
                num_trials = length(cfg.trl);%number of trials
                for t = 1:num_trials
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
                        
                        LFPs = data(LFPchannels(chan)).values{t};
                        for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                            fixt = fixationtimes(1,f);%start of fixation
                            fixdur = fixationtimes(2,f)-fixationtimes(1,f);
                            prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                            if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                                continue;
                            end
                            sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                            if sacamp < min_sac_amp %prior saccade is too small so ignore
                                continue
                            end
                                                        
                            %get firing rate locked to fixation
                            this_fixLFP = LFPs(fixt-twin1+1:fixt+fix_dur_to_warp_to);
                            if fixdur < fix_dur_to_warp_to && cutoff_short_fixations
                                remaining = fix_dur_to_warp_to-fixdur;
                                this_fixLFP(end-remaining+1:end) = NaN;
                            end
                            fix_LFP(fix_ind,:) = this_fixLFP;
                            fix_LFP_RMS(fix_ind) = sqrt(nansum((this_fixLFP(twin1-t_start+1:end)).^2));

                                
                            
                            if fixdur == fix_dur_to_warp_to
                                fix_LFP_resampled(fix_ind,:) = this_fixLFP;
                            else
                                %resampling
                                temp = NaN(1,twin1+fix_dur_to_warp_to);
                                temp(1:twin1) = this_fixLFP(1:twin1);
                                temp2 = LFPs(fixt+1:fixt+fixdur);
                                resampled = round(linspace(1,fixdur,fix_dur_to_warp_to));
                                temp2 = temp2(resampled);
                                temp(twin1+1:end) = temp2;
                                fix_LFP_resampled(fix_ind,:) = temp;           
                            end
                            fix_LFP_RMS_resampled(fix_ind) = sqrt(nansum((temp(twin1-t_start+1:end)).^2));

                            all_LFP_fixdur(fix_ind) = fixdur;
                            fix_ind = fix_ind+1;
                        end
                    end
                end
                
                %---Remove Excess NaNs---%
                fix_LFP = laundry(fix_LFP);
                fix_LFP_resampled = laundry(fix_LFP_resampled);
                all_LFP_fixdur = laundry(all_LFP_fixdur);
                
                if plot_individual_cells
                    [~,si] = sort(all_LFP_fixdur);
                    fix_LFP = fix_LFP(si,:);
                    fix_LFP_resampled(si,:);
                    sm_fix_LFP = imgaussfilt(fix_LFP,8);
                    sm_LFP_warped = imgaussfilt(fix_LFP_resampled, 8);
                    
                    figure
                    subplot(2,2,1)
                    imagesc(-twin1:fix_dur_to_warp_to-1,1:size(fix_LFP,1),sm_fix_LFP)
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Fixation #')
                    title('Raw-Sorted by Fixation Duration')
                    box off
                    
                    subplot(2,2,3)
                    imagesc(-twin1:fix_dur_to_warp_to-1,1:size(fix_LFP_resampled,1),sm_LFP_warped)
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Fixation #')
                    title('Warped-Sorted by Fixation Duration')
                    box off
                    
                    subplot(2,2,2)
                    plot(-twin1:fix_dur_to_warp_to-1,nanmean(fix_LFP),'k')
                    hold on
                    plot(-twin1:fix_dur_to_warp_to-1,nanmean(fix_LFP_resampled),'b')
                    plot([-twin1 fix_dur_to_warp_to],[0 0],'k--')
                    yl = ylim;
                    plot([0 0],[yl(1) yl(2)],'k--')
                    hold off
                    xlim([-twin1 fix_dur_to_warp_to])
                    legend('Raw','Resampled')
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Fixation #')
                    box off
                    
                    yraw = nanmean(fix_LFP(:,twin1-50:end));
                    yresampled = nanmean(fix_LFP_resampled(:,twin1-50:end));
                    
                    lambda_x = yraw-min(yraw);%lambda_x must be positive
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    raw_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    lambda_x = yresampled-min(yresampled);%lambda_x must be positive
                    lambda = nansum(nansum(lambda_x.*p_x));
                    plogp = lambda_x.*log2(lambda_x/lambda);
                    resampled_skaggs = nansum(nansum(plogp.*p_x)); %bits/second
                    
                    fixation_LFP_skagg_info_raw_resampled(1,LFP_ind) = raw_skaggs;
                    fixation_LFP_skagg_info_raw_resampled(2,LFP_ind) = resampled_skaggs;
                    
                    fixation_LFP_amplitude_resampled(1,LFP_ind) = nanmean(fix_LFP_RMS);
                    fixation_LFP_amplitude_resampled(2,LFP_ind) = nanmean(fix_LFP_RMS_resampled);
                    
                    p_change_LFP = 100*(resampled_skaggs-raw_skaggs)/raw_skaggs;
                    
                    title(['All Fixations: ' num2str(p_change_LFP,3) '% Change'])
                    
                    subtitle([task_file(1:end-11) ' LFP Channel ' num2str(chan)])
                    save_and_close_fig(figure_dir2,[task_file(1:end-11) '_LFP_Channel' num2str(chan) '_LFP_fixation_time_warping']);
                end
                %%
                LFP_ind = LFP_ind +1;
            end
        end
    end
end
set(0,'DefaultFigureVisible','ON');
%% Change in in Fixation-aligned Temporal Information with Resampling

%---Out2In Fixations
[~,p_fix_out2in] = ttest(fixation_skagg_info_raw_resampled(1,:),fixation_skagg_info_raw_resampled(2,:));
change_fix_out2in = 100*(fixation_skagg_info_raw_resampled(2,:)-fixation_skagg_info_raw_resampled(1,:))./fixation_skagg_info_raw_resampled(1,:);
median_change_fix_out2in = nanmedian(change_fix_out2in);
change_fix_out2in(change_fix_out2in < -300) = -300;%for ploting make it nicer looking
change_fix_out2in(change_fix_out2in > 300) = 300;%for ploting make it nicer looking

[~,p_fr_out2in] = ttest(fixation_average_firing_rate_resampled(1,:),fixation_average_firing_rate_resampled(2,:));
fr_change_fix_out2in = 100*(fixation_average_firing_rate_resampled(2,:)-fixation_average_firing_rate_resampled(1,:))./fixation_average_firing_rate_resampled(1,:);
median_fr_change_fix_out2in = nanmedian(fr_change_fix_out2in);
fr_change_fix_out2in(fr_change_fix_out2in < -200) = 200;
fr_change_fix_out2in(fr_change_fix_out2in > 200) = 200;

vals1 = change_fix_out2in;
vals1(isnan(vals1)) = [];
vals2 = fr_change_fix_out2in;
vals2(isnan(vals2)) = [];
corr_change_fix_out2in = corrcoef(vals1,vals2);

figure
subplot(2,3,1)
hist(change_fix_out2in,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['Out2in Only: median change = ' num2str(median_change_fix_out2in,3) '%, (p = ' num2str(p_fix_out2in,3) ')'])

subplot(2,3,2)
hist(fr_change_fix_out2in,25)
box off
xlabel('% Change in Firing Rate')
ylabel('View Cell Count')
title(['Out2in Only: median change = ' num2str(median_fr_change_fix_out2in,3) '%, (p = ' num2str(p_fr_out2in,3) ')'])

subplot(2,3,3)
plot(change_fix_out2in,fr_change_fix_out2in,'.k')
xlabel('% Change in Skagg')
ylabel('% Change in Firing Rate')
title(['Out2In: Corr Change in Firing Rate & Skagg, r = ' num2str(corr_change_fix_out2in(2),3)])
box off

%---All Fixations---%
[~,p_fix_all] = ttest(fixation_skagg_info_raw_resampled(3,:),fixation_skagg_info_raw_resampled(4,:));
change_fix_all = 100*(fixation_skagg_info_raw_resampled(4,:)-fixation_skagg_info_raw_resampled(3,:))./fixation_skagg_info_raw_resampled(3,:);
mean_change_fix_all = nanmedian(change_fix_all);
change_fix_all(change_fix_all < -300) = -300;%for ploting make it nicer looking
change_fix_all(change_fix_all > 300) = 300;%for ploting make it nicer looking

[~,p_fr_all] = ttest(fixation_average_firing_rate_resampled(3,:),fixation_average_firing_rate_resampled(4,:));
fr_change_fix_all = 100*(fixation_average_firing_rate_resampled(4,:)-fixation_average_firing_rate_resampled(3,:))./fixation_average_firing_rate_resampled(3,:);
median_fr_change_fix_all = nanmedian(fr_change_fix_all);
fr_change_fix_all(fr_change_fix_all < -200) = 200;
fr_change_fix_all(fr_change_fix_all > 200) = 200;

vals1 = change_fix_all;
vals1(isnan(vals1)) = [];
vals2 = fr_change_fix_all;
vals2(isnan(vals2)) = [];
corr_change_fix_all = corrcoef(vals1,vals2);

subplot(2,3,4)
hist(change_fix_all,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(mean_change_fix_all,3) '%, (p = ' num2str(p_fix_all,3) ')'])

subplot(2,3,5)
hist(fr_change_fix_all,25)
box off
xlabel('% Change in Firing Rate')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(median_fr_change_fix_all,3) '%, (p = ' num2str(p_fr_all,3) ')'])

subplot(2,3,6)
plot(change_fix_all,fr_change_fix_all,'.k')
xlabel('% Change in Skagg')
ylabel('% Change in Firing Rate')
title(['All Fixations: Corr Change in Firing Rate & Skagg, r = ' num2str(corr_change_fix_all(2),3)])
box off

%% Change in in Fixation-aligned Temporal Information with Scaling

%---Out2In Fixations
[~,p_scaled_fix_out2in] = ttest(fixation_skagg_info_raw_scaled(1,:),fixation_skagg_info_raw_scaled(2,:));
change_scaled_fix_out2in = 100*(fixation_skagg_info_raw_scaled(2,:)-fixation_skagg_info_raw_scaled(1,:))./fixation_skagg_info_raw_scaled(1,:);
median_change_scaled_fix_out2in = nanmedian(change_scaled_fix_out2in);
change_scaled_fix_out2in(change_scaled_fix_out2in < -300) = -300;%for ploting make it nicer looking
change_scaled_fix_out2in(change_scaled_fix_out2in > 300) = 300;%for ploting make it nicer looking

[~,p_scaled_fr_out2in] = ttest(fixation_average_firing_rate_scaled(1,:),fixation_average_firing_rate_scaled(2,:));
fr_change_scaled_fix_out2in = 100*(fixation_average_firing_rate_scaled(2,:)-fixation_average_firing_rate_scaled(1,:))./fixation_average_firing_rate_scaled(1,:);
median_fr_change_scaled_fix_out2in = nanmedian(fr_change_scaled_fix_out2in);
fr_change_scaled_fix_out2in(fr_change_scaled_fix_out2in < -200) = 200;
fr_change_scaled_fix_out2in(fr_change_scaled_fix_out2in > 200) = 200;

vals1 = change_scaled_fix_out2in;
vals1(isnan(vals1)) = [];
vals2 = fr_change_scaled_fix_out2in;
vals2(isnan(vals2)) = [];
corr_change_scaled_fix_out2in = corrcoef(vals1,vals2);

figure
subplot(2,3,1)
hist(change_scaled_fix_out2in,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['Out2in Only: median change = ' num2str(median_change_scaled_fix_out2in,3) '%, (p = ' num2str(p_scaled_fix_out2in,3) ')'])

subplot(2,3,2)
hist(fr_change_scaled_fix_out2in,25)
box off
xlabel('% Change in Firing Rate')
ylabel('View Cell Count')
title(['Out2in Only: median change = ' num2str(median_fr_change_scaled_fix_out2in,3) '%, (p = ' num2str(p_scaled_fr_out2in,3) ')'])

subplot(2,3,3)
plot(change_scaled_fix_out2in,fr_change_scaled_fix_out2in,'.k')
xlabel('% Change in Skagg')
ylabel('% Change in Firing Rate')
title(['Out2In: Corr Change in Firing Rate & Skagg, r = ' num2str(corr_change_scaled_fix_out2in(2),3)])
box off

%---All Fixations---%
[~,p_scaled_fix_all] = ttest(fixation_skagg_info_raw_scaled(3,:),fixation_skagg_info_raw_scaled(4,:));
change_scaled_fix_all = 100*(fixation_skagg_info_raw_scaled(4,:)-fixation_skagg_info_raw_scaled(3,:))./fixation_skagg_info_raw_scaled(3,:);
mean_change_scaled_fix_all = nanmedian(change_scaled_fix_all);
change_scaled_fix_all(change_scaled_fix_all < -300) = -300;%for ploting make it nicer looking
change_scaled_fix_all(change_scaled_fix_all > 300) = 300;%for ploting make it nicer looking

[~,p_scaled_fr_all] = ttest(fixation_average_firing_rate_scaled(3,:),fixation_average_firing_rate_scaled(4,:));
fr_change_scaled_fix_all = 100*(fixation_average_firing_rate_scaled(4,:)-fixation_average_firing_rate_scaled(3,:))./fixation_average_firing_rate_scaled(3,:);
median_fr_change_scaled_fix_all = nanmedian(fr_change_scaled_fix_all);
fr_change_scaled_fix_all(fr_change_scaled_fix_all < -200) = 200;
fr_change_scaled_fix_all(fr_change_scaled_fix_all > 200) = 200;

vals1 = change_scaled_fix_all;
vals1(isnan(vals1)) = [];
vals2 = fr_change_scaled_fix_all;
vals2(isnan(vals2)) = [];
corr_change_scaled_fix_all = corrcoef(vals1,vals2);

subplot(2,3,4)
hist(change_scaled_fix_all,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(mean_change_scaled_fix_all,3) '%, (p = ' num2str(p_scaled_fix_all,3) ')'])

subplot(2,3,5)
hist(fr_change_scaled_fix_all,25)
box off
xlabel('% Change in Firing Rate')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(median_fr_change_scaled_fix_all,3) '%, (p = ' num2str(p_scaled_fr_all,3) ')'])

subplot(2,3,6)
plot(change_scaled_fix_all,fr_change_scaled_fix_all,'.k')
xlabel('% Change in Skagg')
ylabel('% Change in Firing Rate')
title(['All Fixations: Corr Change in Firing Rate & Skagg, r = ' num2str(corr_change_scaled_fix_all(2),3)])
box off

%% LFP Fixation Aligned Responses

%---All Fixations---%
[~,p_fix_LFP_] = ttest(fixation_LFP_skagg_info_raw_resampled(1,:),fixation_LFP_skagg_info_raw_resampled(2,:));
change_fix_LFP_ = 100*(fixation_LFP_skagg_info_raw_resampled(2,:)-fixation_LFP_skagg_info_raw_resampled(1,:))./fixation_LFP_skagg_info_raw_resampled(1,:);
mean_change_fix_LFP_ = nanmedian(change_fix_LFP_);

[~,p_LFP_all] = ttest(fixation_LFP_amplitude_resampled(1,:),fixation_LFP_amplitude_resampled(2,:));
LFP_change_fix_LFP_ = 100*(fixation_LFP_amplitude_resampled(2,:)-fixation_LFP_amplitude_resampled(1,:))./fixation_LFP_amplitude_resampled(1,:);
mean_LFP_change_fix_LFP_ = nanmedian(LFP_change_fix_LFP_);

vals1 = change_fix_LFP_;
vals1(isnan(vals1)) = [];
vals2 = LFP_change_fix_LFP_;
vals2(isnan(vals2)) = [];
corr_change_fix_LFP_ = corrcoef(vals1,vals2);

figure
subplot(2,2,1)
hist(change_fix_LFP_,25)
box off
xlabel('% Change in Skagg Info')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(mean_change_fix_LFP_,3) '%, (p = ' num2str(p_fix_LFP_,3) ')'])

subplot(2,2,2)
hist(LFP_change_fix_LFP_,25)
box off
xlabel('% Change in LFP Amplitude (RMS)')
ylabel('View Cell Count')
title(['All Fixations: median change = ' num2str(mean_LFP_change_fix_LFP_,3) '%, (p = ' num2str(p_LFP_all,3) ')'])

subplot(2,2,3)
plot(change_fix_LFP_,LFP_change_fix_LFP_,'.k')
xlabel('% Change in Skagg')
ylabel('% Change in Firing Rate')
title(['All Fixations: Corr Change in Firing Rate & Skagg, r = ' num2str(corr_change_fix_LFP_(2),3)])
box off        