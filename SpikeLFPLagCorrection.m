function SpikeLFPLagCorrection(data_dir,figure_dir,session_data)
%Written by Seth Konig 2/18/19 to correct for a data mis-alignment/lag between
%Spike times and analog data. Analag data including eye data and LFPs
%appears to be up to 70 ms behind the spike times. The event codes appear
%to the best of our knowledge synched well with the spike times. The eye
%data is too hard to use for alignment so use spike-LFP alignment since
%very consitent and easy to measure

lag_threshold4Correction = 20;%ms anything less does not need correction
median_trough_time = 5;%ms what is the expected measure for the trough when no apparent lag
median_peak_time = 9;%ms what is the expected measure for the peak when no apparent lag so subtract for correction

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

%load expected image evoked LFP when no apparent lag
load('C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\ImageEvokedLFPNoLag.mat');

lag_correction_data = [];
lag_correction_data.lag_threshold4Correction = lag_threshold4Correction;
lag_correction_data.median_trough_time = median_trough_time;
lag_correction_data.median_peak_time = median_peak_time;
lag_correction_data.correction_applied = 0;

task = 'ListSQ';
Fs = 1000; %Hz sampling frequency
twin = 250;
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer


figure_dir = [figure_dir 'SpikeLFPLagCorrection\'];
[bhigh,ahigh] = butter(8,20/(Fs/2),'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read in task file data
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,...
    sorting_quality,waveform_count,lfp_quality,comments] = get_task_data(session_data,task);
if isempty(task_file)
    return
end

%load task file data
load([data_dir task_file(1:end-11) '-preprocessed.mat']);

%grab unit data
[~,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    save([data_dir task_file(1:end-11) '-preprocessed.mat'],'lag_correction_data','-append')
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

num_trials = length(cfg.trl); %number of trials
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    save([data_dir task_file(1:end-11) '-preprocessed.mat'],'lag_correction_data','-append')
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return;
end

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

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
%LFPchannels(bad_channels) = NaN;
% if all(isnan(LFPchannels))
%     disp([task_file(1:8) ': all LFP channels are bad ... be cautious'])
% end

%---Get Spike-Triggered Average LFP to determine Correction---%
pre_corrected_STA_LFP = cell(1,num_units);
for unit = 1:num_units
    unit_channel = str2double(unit_stats{1,unit}(6));
    if unit_channel > length(LFPchannels)
        continue
    end
    pre_corrected_STA_LFP{unit} = NaN(length(waveforms{2,unit}),twin*2+1);
    spike_ind = 1;
    
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            trial_start = cfg.trl(t).alltim(1);
            trial_end = cfg.trl(t).alltim(end)-trial_start;
            
            spikes = find(data(unit).values{t}); %spike trains for this trial
            spikes(spikes <= twin) = [];
            spikes(spikes >= trial_end-twin) = [];
            
            LFPs = data(LFPchannels(unit_channel)).values{t};
            spikes(spikes > length(LFPs)-twin) = [];
            for spk = 1:length(spikes);
                pre_corrected_STA_LFP{unit}(spike_ind,:) = LFPs(spikes(spk)-twin:spikes(spk)+twin);
                spike_ind = spike_ind+1;
            end
        end
    end
end
pre_corrected_STA_LFP = laundry(pre_corrected_STA_LFP);

%---Get Image Evoked LFP for extra verification---%
img_ind = 1;
pre_corrected_image_ta_LFP = cell(1,length(LFPchannels));
for chan = 1:length(LFPchannels)
    pre_corrected_image_ta_LFP{chan} = NaN(196,twin*2+1);
end
for t = 1:num_trials
    if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
        imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
        for chan = 1:length(LFPchannels)
            LFPs = data(LFPchannels(chan)).values{t};
            pre_corrected_image_ta_LFP{chan}(img_ind,:) = LFPs(imgon-twin:imgon+twin);
        end
        img_ind = img_ind+1;
    end
end
pre_corrected_image_ta_LFP = laundry(pre_corrected_image_ta_LFP);

%---Get Saccade Triggered Average LFP to Make sure I don't mess with Eye-LFP relationship---%
%yes I purposely go through each unit even though it's really by
%channel...just want to look at things a little differently and spike-eye
%relationships are super important to observe!
pre_sacTA_LFP = cell(1,num_units);
for unit = 1:num_units
    pre_sacTA_LFP{unit} = zeros(1,2*twin+1);
    fix_count = 0;
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
                invalid= find(fixationtimes(1,:) > imgoff-twin);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %saccade started before image turned on
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(saccadetimes(1,:) > imgoff-twin);
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
                
                if sum(~isnan(LFPchannels)) == 1 %chose the only available channel
                    LFPs = data(LFPchannels(~isnan(LFPchannels))).values{t};
                else
                    unit_channel = str2double(unit_stats{1,unit}(6));
                    if unit_channel > length(LFPchannels)
                        first_non_nan = find(~isnan(LFPchannels));
                        first_non_nan = first_non_nan(1);
                        LFPs = data(LFPchannels(first_non_nan)).values{t};
                    else
                        if ~isnan(LFPchannels(unit_channel))
                            LFPs = data(LFPchannels(unit_channel)).values{t};
                        else
                            first_non_nan = find(~isnan(LFPchannels));
                            first_non_nan = first_non_nan(1);
                            LFPs = data(LFPchannels(first_non_nan)).values{t};
                        end
                    end
                end
                
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
                    fixt = fixationtimes(1,f);%start of fixation
                    sact = saccadetimes(1,prior_sac);
                    
                    pre_sacTA_LFP{unit} = pre_sacTA_LFP{unit}+LFPs(sact-twin:sact+twin);
                    fix_count = fix_count+1;
                end
            end
        end
    end
    pre_sacTA_LFP{unit} = pre_sacTA_LFP{unit}/fix_count;
end


%%
%---High Pass STA-LFP and Image Evoked Response---%
avg_pre_corrected_STA_LFP = cell(1,num_units);
for unit = 1:num_units
    if ~isempty(pre_corrected_STA_LFP{unit})
        avg_STA_LFP = nanmean(pre_corrected_STA_LFP{unit});
        avg_pre_corrected_STA_LFP{unit} = filtfilt(bhigh,ahigh,avg_STA_LFP);
    end
end

avg_pre_corrected_image_ta_LFP = NaN(length(LFPchannels),2*twin+1);
for chan = 1:length(LFPchannels)
    avg_STA_LFP = mean(pre_corrected_image_ta_LFP{chan});
    avg_pre_corrected_image_ta_LFP(chan,:) = filtfilt(bhigh,ahigh,avg_STA_LFP);
end

%---Determine Lag from STA-LFP---%
all_peaks_all_troughs = NaN(4,num_units);
for unit = 1:num_units
    if ~isempty(avg_pre_corrected_STA_LFP{unit})
        rms = sqrt(mean(avg_pre_corrected_STA_LFP{unit}.^2));
        [PKS,LOCS] = findpeaks(avg_pre_corrected_STA_LFP{unit},'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
        [trofs,trLOCS] = findpeaks(-avg_pre_corrected_STA_LFP{unit},'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
        
        if length(LOCS) == 1 %only want well detected ones
            all_peaks_all_troughs(1,unit) = LOCS;
            all_peaks_all_troughs(2,unit) = PKS;
        end
        if length(trLOCS) == 1 %only want well detected ones
            all_peaks_all_troughs(3,unit) = trLOCS;
            all_peaks_all_troughs(4,unit) = trofs;
        end
    end
end
%% Check Variance and mean troughs and peaks
var_trough = nanstd(all_peaks_all_troughs(3,:)-twin);
mean_trough = nanmean(all_peaks_all_troughs(3,:)-twin);
if var_trough > 2
    disp([task_file(1:end-11) 'Lots of Variance in when trough was detected']);
    trough_correction_to_apply = 0;
else
    if mean_trough > 20
        trough_correction_to_apply = round(mean_trough-median_trough_time);
    else
        trough_correction_to_apply = 0;
    end
end
var_peak = nanstd(all_peaks_all_troughs(1,:)-twin);
mean_peak = nanmean(all_peaks_all_troughs(1,:)-twin);
if var_peak > 2
    disp([task_file(1:end-11) 'Lots of Variance in when peak was detected'])
else
    if mean_peak > 20
        peak_correction_to_apply = round(mean_peak-median_peak_time);
    else
        peak_correction_to_apply = 0;
    end
end

if peak_correction_to_apply > 0 || trough_correction_to_apply > 0
    if abs(trough_correction_to_apply-peak_correction_to_apply) < 4 %very similar
        correction_to_apply = trough_correction_to_apply;
    else
        disp('Correction to Apply for Peak and trough are very diffrent!')
        correction_to_apply = trough_correction_to_apply;
    end
else
    correction_to_apply = 0;
end
%% Apply correction to LFPs

if correction_to_apply > 0
    %---Shift EYe Data by correction and apply clusterfix---%
    eyechans = find_desired_channels(cfg,'eye');
    new_eye_data = cell(1,length(eyechans));
    for chan = 1:length(eyechans)
        new_eye_data{chan}.values = cell(1,num_trials);
        for t = 1:num_trials
            eyedt = data(eyechans(chan)).values{t};
            if t ~= num_trials
                eyedt2 = data(eyechans(chan)).values{t+1}(1:correction_to_apply);
            else
                eyedt2 = [];
            end
            new_eye_data{chan}.values{t} = [eyedt(correction_to_apply+1:end) eyedt2];
        end
    end
    
    %Running Cluster Fix
    new_fixationstats = cell(1,num_trials);
    for t = 1:num_trials;
        if ~isempty(new_eye_data{1}.values{t})
            new_fixationstats{t} = ClusterFix_Plexon([new_eye_data{1}.values{t};new_eye_data{2}.values{t}]);
        end
    end
    
    %---Shift LFPs by correction---%
    new_LFP_data = cell(1,length(LFPchannels));
    for chan = 1:length(LFPchannels)
        new_LFP_data{chan}.values = cell(1,num_trials);
        for t = 1:num_trials
            LFPs = data(LFPchannels(chan)).values{t};
            if t ~= num_trials
                LFPs2 = data(LFPchannels(chan)).values{t+1}(1:correction_to_apply);
            else
                LFPs2 = [];
            end
            new_LFP_data{chan}.values{t} = [LFPs(correction_to_apply+1:end) LFPs2];
        end
    end
    
    %---Get Spike-Triggered Average LFP to determine Correction---%
    post_corrected_STA_LFP = cell(1,num_units);
    for unit = 1:num_units
        unit_channel = str2double(unit_stats{1,unit}(6));
        if unit_channel > length(LFPchannels)
            first_non_nan = find(~isnan(LFPchannels));
            unit_channel = first_non_nan(1);
        end
        post_corrected_STA_LFP{unit} = NaN(length(waveforms{2,unit}),twin*2+1);
        spike_ind = 1;
        
        for t = 1:num_trials
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                trial_start = cfg.trl(t).alltim(1);
                trial_end = cfg.trl(t).alltim(end)-trial_start;
                
                spikes = find(data(unit).values{t}); %spike trains for this trial
                spikes(spikes <= twin) = [];
                spikes(spikes >= trial_end-twin) = [];
                
                LFPs = new_LFP_data{unit_channel}.values{t};
                if t ~= num_trials
                    if length(LFPs) < max(spikes)
                        error('Should have LFP data when spikes occured')
                    end
                else
                    if length(LFPs) < max(spikes)
                        disp('LFP shorted than max spike count')
                    end
                    spikes(spikes > length(LFPs)-twin) = [];
                end
                for spk = 1:length(spikes);
                    post_corrected_STA_LFP{unit}(spike_ind,:) = LFPs(spikes(spk)-twin:spikes(spk)+twin);
                    spike_ind = spike_ind+1;
                end
            end
        end
    end
    post_corrected_STA_LFP = laundry(post_corrected_STA_LFP);
    
    %---Get Image Evoked LFP for extra verification---%
    img_ind = 1;
    post_corrected_image_ta_LFP = cell(1,length(LFPchannels));
    for chan = 1:length(LFPchannels)
        post_corrected_image_ta_LFP{chan} = NaN(196,twin*2+1);
    end
    for t = 1:num_trials
        if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
            for chan = 1:length(LFPchannels)
                LFPs = new_LFP_data{unit_channel}.values{t};
                post_corrected_image_ta_LFP{chan}(img_ind,:) = LFPs(imgon-twin:imgon+twin);
            end
            img_ind = img_ind+1;
        end
    end
    post_corrected_image_ta_LFP = laundry(post_corrected_image_ta_LFP);
    
    %---Average and High Pass---%
    avg_post_corrected_STA_LFP = cell(1,num_units);
    for unit = 1:num_units
        if ~isempty(pre_corrected_STA_LFP{unit})
            avg_STA_LFP = nanmean(post_corrected_STA_LFP{unit});
            avg_post_corrected_STA_LFP{unit} = filtfilt(bhigh,ahigh,avg_STA_LFP);
        end
    end
    
    avg_post_corrected_image_ta_LFP = NaN(length(LFPchannels),2*twin+1);
    for chan = 1:length(LFPchannels)
        avg_STA_LFP = mean(post_corrected_image_ta_LFP{chan});
        avg_post_corrected_image_ta_LFP(chan,:) = filtfilt(bhigh,ahigh,avg_STA_LFP);
    end
    
    %---Get Saccade Triggered Average LFP to Make sure I don't mess with Eye-LFP relationship---%
    %yes I purposely go through each unit even though it's really by
    %channel...just want to look at things a little differently and spike-eye
    %relationships are super important to observe!
    post_sacTA_LFP = cell(1,num_units);
    for unit = 1:num_units
        post_sacTA_LFP{unit} = zeros(1,2*twin+1);
        fix_count = 0;
        for t = 1:num_trials
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                    imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                    imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                    
                    %---image info---%
                    img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                    
                    % if monkey isn't paying attention and looked away image postsentation
                    % is now longer than imgdur (because of cumulative looking time)
                    % so data isn't probably worth much plus have to cut off somewhere
                    if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                        imgoff = imgon+1.5*imgdur-1;
                    end
                    imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                    
                    fixationtimes = new_fixationstats{t}.fixationtimes; %fixation start and end times
                    saccadetimes = new_fixationstats{t}.saccadetimes; %saccade start and end times
                    fixations = round(new_fixationstats{t}.fixations); %mean fixation location
                    xy = new_fixationstats{t}.XY; %xy eye trace
                    
                    %find fiations and saccades that did not occur during the image period;
                    %should also take care of the 1st fixation on the crosshair
                    
                    %fixation started before image turned on
                    invalid= find(fixationtimes(1,:) < imgon);
                    fixationtimes(:,invalid) = [];
                    fixations(:,invalid) = [];
                    
                    %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                    invalid= find(fixationtimes(1,:) > imgoff-twin);
                    fixationtimes(:,invalid) = [];
                    fixations(:,invalid) = [];
                    
                    %saccade started before image turned on
                    invalid= find(saccadetimes(1,:) < imgon);
                    saccadetimes(:,invalid) = [];
                    
                    %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                    invalid= find(saccadetimes(1,:) > imgoff-twin);
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
                    
                    if sum(~isnan(LFPchannels)) == 1 %chose the only available channel
                        LFPs = data(LFPchannels(~isnan(LFPchannels))).values{t};
                    else
                        unit_channel = str2double(unit_stats{1,unit}(6));
                        if unit_channel > length(LFPchannels)
                            first_non_nan = find(~isnan(LFPchannels));
                            first_non_nan = first_non_nan(1);
                            LFPs = new_LFP_data{first_non_nan}.values{t};
                        else
                            if ~isnan(LFPchannels(unit_channel))
                                LFPs = new_LFP_data{unit_channel}.values{t};
                            else
                                first_non_nan = find(~isnan(LFPchannels));
                                first_non_nan = first_non_nan(1);
                                LFPs = new_LFP_data{first_non_nan}.values{t};
                            end
                        end
                    end
                    
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
                        fixt = fixationtimes(1,f);%start of fixation
                        sact = saccadetimes(1,prior_sac);
                        
                        post_sacTA_LFP{unit} = post_sacTA_LFP{unit}+LFPs(sact-twin:sact+twin);
                        fix_count = fix_count+1;
                    end
                end
            end
        end
        post_sacTA_LFP{unit} = post_sacTA_LFP{unit}/fix_count;
    end
else
    avg_post_corrected_STA_LFP = [];
    avg_post_corrected_image_ta_LFP = [];
    post_sacTA_LFP = [];
end
%%

figure
subplot(2,2,1)
unit = find(~cellfun(@isempty,pre_corrected_STA_LFP));%just plot 1st unit available for example
unit = unit(1);
pct5 = prctile(pre_corrected_STA_LFP{unit}(:),1);
pct95 = prctile(pre_corrected_STA_LFP{unit}(:),99);
imagesc(-twin:twin,1:size(pre_corrected_STA_LFP{unit},1),pre_corrected_STA_LFP{unit})
caxis([pct5 pct95])
xlim([-twin/4 twin/2])
xlabel('Time from Spike (ms)')
ylabel('Spike #')
title(['Pre-Corrected Raw spike-aligned LFP for unit ' num2str(unit)])
box off

subplot(2,2,2)
hold on
plot([-twin twin],[0 0],'k')
p = [];
pc = [];
pind = 1;
for unit = 1:num_units
    if ~isempty(avg_pre_corrected_STA_LFP{unit})
        p(pind) = plot(-twin:twin,avg_pre_corrected_STA_LFP{unit}+unit,'b');
        if ~isempty(avg_post_corrected_STA_LFP)
            pc(pind) = plot(-twin:twin,avg_post_corrected_STA_LFP{unit}+num_units+unit,'r');
        end
        pind = pind+1;
    end
end
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-twin/4 twin/2])
xlabel('Time from Spike (ms)')
ylabel('Average LFP (uV)')
if ~isempty(pc)
    legend([p(1) pc(1)],'Pre-Corrected','Post-Corrected')
else
    legend(p(1),'Pre-Corrected')
end
title(['STA-LFP: Trough @ ' num2str(mean_trough,2) ' ms, Peak @ ' num2str(mean_peak,2) ' ms'])

subplot(2,2,3)
hold on
p2 = [];
p2c = [];
pind = 1;
plot([-twin twin],[0 0],'k')
for unit = 1:num_units
    if ~isempty(pre_sacTA_LFP{unit})
        p2(pind) = plot(-twin:twin,pre_sacTA_LFP{unit}+unit,'b');
        if ~isempty(post_sacTA_LFP)
            p2c(pind) = plot(-twin:twin,post_sacTA_LFP{unit}+unit+num_units,'r');
        end
        pind = pind+1;
    end
end
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-twin twin])
xlabel('Time from Spike (ms)')
ylabel('Average LFP (uV)')
if ~isempty(p2c)
    legend([p2(1) p2c(1)],'Pre-Corrected','Post-Corrected')
else
    legend(p2(1),'Pre-Corrected')
end
title('Saccade TA-LFP-Should Not change with Lag Correction')

subplot(2,2,4)
p3 = [];
pc3 = [];
pind = 1;
hold on
plot([-twin twin],[0 0],'k')
for chan = 1:length(LFPchannels)
    p3(pind) = plot(-twin:twin,avg_pre_corrected_image_ta_LFP(chan,:)-mean(avg_pre_corrected_image_ta_LFP(chan,:))+chan,'b');
    if ~isempty(avg_post_corrected_image_ta_LFP)
        pc3(pind) = plot(-twin:twin,avg_post_corrected_image_ta_LFP(chan,:)-mean(avg_post_corrected_image_ta_LFP(chan,:))+chan,'r');
    end
    pind = pind+1;
end
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
avgIE = avgImageEvokedLFP-mean(avgImageEvokedLFP);
avgIE = (avgIE/sqrt(sum(avgIE.^2)))*sqrt(sum(mean(avg_pre_corrected_image_ta_LFP).^2)); %RMS mean
p4 = plot(-twin:twin,avgIE,'g');
hold off
xlim([-twin/4 twin])
xlabel('Time from Image Onset (ms)')
ylabel('[Scaled-zeroed] Average LFP (uV)')
if ~isempty(pc3)
    legend([p3(1) pc3(1) p4(1)],'Pre-Corrected','Post-Corrected','"Good" Pop. Avg.')
else
    legend([p3(1) p4(1)],'Pre-Corrected','"Good" Pop. Avg.')
end
title('Image Evoked LFP-Should Change with Lag Correction')

subtitle([task_file(1:end-11) ', Applying correction of ' num2str(correction_to_apply) ' ms'])
%%
save_and_close_fig(figure_dir,[task_file(1:end-11) '_LagCorrection'])

%%
%%---Save Data---%
lag_correction_data.correction_applied = correction_to_apply;
lag_correction_data.observed_peak = mean_peak;
lag_correction_data.observed_peak_var = var_peak;
lag_correction_data.observed_trough = mean_trough;
lag_correction_data.observed_trough_var = var_trough;
lag_correction_data.observed_STA_LFP = pre_corrected_STA_LFP;
lag_correction_data.observed_sacTA_LFP = pre_sacTA_LFP;
lag_correction_data.observed_image_ta_LFP = avg_pre_corrected_image_ta_LFP;
%will be blank if no correction
lag_correction_data.corrected_STA_LFP = avg_post_corrected_STA_LFP;
lag_correction_data.corrected_sacTA_LFP = post_sacTA_LFP;
lag_correction_data.corrected_image_ta_LFP = avg_post_corrected_image_ta_LFP;

if correction_to_apply > 0
    %re-assigned shifted eye data, fixationstats, and LFPs
    fixationstats = new_fixationstats;
    for chan = 1:length(eyechans)
        data(eyechans(chan)).values = new_eye_data{chan}.values;
    end
    for chan = 1:length(LFPchannels)
        data(LFPchannels(chan)).values = new_LFP_data{chan}.values;
    end
    
    save([data_dir task_file(1:end-11) '-preprocessed.mat'],'lag_correction_data','fixationstats','data','-append')
else
    save([data_dir task_file(1:end-11) '-preprocessed.mat'],'lag_correction_data','-append')
end

end