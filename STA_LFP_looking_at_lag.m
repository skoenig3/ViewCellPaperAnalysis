% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig February 2019
% Code does the following
% 1) Summarizes spatial correlations for place cell and non-place cells
% 2) Calculates fixation aligned firing rate curves for place cells
% 3) Calculates eye coverage and place field coverage for place cells
% 4) Tracks AP location, unit counts, and which monkey (not currently used)
% 5) Contextual differences between list and sequence task
% 6) Copies relevant figures for place cells to summary directory

%Code rechecked by SDK on 1/5/2017
%Small bug fixed on 2/6/18 on designating inside and outside sequence
%trials SDK

clar %clear,clc

task = 'ListSQ';
Fs = 1000; %Hz sampling frequency
twin = 250;
twin2 = 25;
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial

[bhigh,ahigh] = butter(8,20/(Fs/2),'high');

cell_ind = 1;

%Spike time comparison variables
spikes_not_found = NaN(1,350);
offset = NaN(2,350);
first_last_spike = NaN(2,350); %includes unsorted
analog_start_end_ts = NaN(2,350);


%STA-LFP variables
recording_id = NaN(1,350);
recording_name = cell(1,350);
min_max_lag = NaN(2,350);
within_sess_variance = NaN(1,350);
view_cell_id = NaN(1,350);
sac_dir_id = NaN(1,350);

%Get average waveform
trough_peak_aligned_STA = zeros(2,2*twin2+1);
tp_unit_count = ones(1,2); 

%Image-TA-LFP
average_image_ta_LFP = NaN(350,2*twin+1);

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 1:2
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
        disp(['Running monkey#' num2str(monk) ', session #' num2str(sess)])
        
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
        load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info');
        load([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'],'mrls');
        
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{sess},task);
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
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
        
        for unit = 1:num_units
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Waveform Comparison across import functions---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_spike_times = [];
            for t = 1:num_trials %use all trials here since looking for spike matches across import types, not spike validity
                spikes = find(data(unit).values{t}); %spike trains for this trial
                spikes = spikes + cfg.trl(t).begsmpind;
                all_spike_times = [all_spike_times spikes];
            end
            
            closest_neighbor = NaN(2,length(all_spike_times));
            for spk = 1:length(all_spike_times)
                closest_neighbor(1,spk) = min(abs(all_spike_times(spk)-waveforms{2,unit}));
                try
                    closest_neighbor(2,spk) = find(all_spike_times(spk) == waveforms{2,unit}+2);
                catch
                    closest_neighbor(2,spk) = NaN;
                end
            end
            
            spikes_not_found(cell_ind) = sum(isnan(closest_neighbor(2,:)));
            offset(1,cell_ind) = mean(closest_neighbor(1,:));
            offset(2,cell_ind) = max(closest_neighbor(1,:));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---STA LFP Analysis from Preporcessed Data---%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            unit_channel = str2double(unit_stats{1,unit}(6));
            if unit_channel > length(LFPchannels)
                continue
            elseif isnan(LFPchannels(unit_channel))
                continue
            end
            
            if exist([data_dir task_file(1:end-4) '.plx']);
                [adfreq, n, ts, fn, ad] = plx_ad([data_dir task_file(1:end-4) '.plx'], unit_channel-1);
                if length(ts) > 1
                    disp(['TS has more than 1 value, ' task_file(1:end-4)])
                    ts = ts(1);
                end
                analog_start_end_ts(1,cell_ind) = round(1000*ts);
                analog_start_end_ts(2,cell_ind) = round(1000*ts)+n;
                
                first_spike = waveforms{2,unit}(1);
                last_spike = waveforms{2,unit}(end);
                
                [n2, npw, ts2, wave] = plx_waves_v([data_dir task_file(1:end-4) '.plx'], unit_channel, 0);%invalid units
                first_invalid = round(ts2(1)*1000);
                last_invalid = round(ts2(end)*1000);
                
                first_last_spike(1,cell_ind) = min(first_spike,first_invalid);
                first_last_spike(2,cell_ind) = max(last_spike,last_invalid);
            end
            
            spike_ind = 1;
            STA_LFP = NaN(length(waveforms{2,unit}),twin*2+1);
            STA_LFPad = NaN(length(waveforms{2,unit}),twin*2+1);
            
            img_ind = 1;
            image_ta_LFP = NaN(196,twin*2+1);
            
            num_trials = length(cfg.trl);%number of trials
            for t = 1:num_trials
                if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                    trial_start = cfg.trl(t).alltim(1);
                    trial_end = cfg.trl(t).alltim(end)-trial_start;
                    
                    spikes = find(data(unit).values{t}); %spike trains for this trial
                    spikes(spikes <= twin) = [];
                    spikes(spikes >= trial_end-twin) = [];
                    
                    LFPs = data(LFPchannels(unit_channel)).values{t};
                    for spk = 1:length(spikes);
                        STA_LFP(spike_ind,:) = LFPs(spikes(spk)-twin:spikes(spk)+twin);
                        
                        %spiked = cfg.trl(t).begsmpind + spikes(spk)-round(ts*1000);
                        %STA_LFPad(spike_ind,:) = ad(spiked-twin:spiked+twin);
                        
                        spike_ind = spike_ind+1;
                    end
                end
                
                if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                    LFPs = data(LFPchannels(unit_channel)).values{t};
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                    imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                    image_ta_LFP(img_ind,:) = LFPs(imgon-twin:imgon+twin);
                    img_ind = img_ind+1;
                end
            end
            %%
            STA_LFP = laundry(STA_LFP);
                        STA_LFPad = laundry(STA_LFPad);
            image_ta_LFP = laundry(image_ta_LFP);
            average_image_ta_LFP(cell_ind,:) = mean(record);
            
            if size(STA_LFP,1) < 100 %% not enough spikes
                continue
            end
            
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                view_cell_id(cell_ind) = 1;
            else
                view_cell_id(cell_ind) = 0;
            end
            if mrls.all_saccades_shuffled_prctile(unit) > 95
                sac_dir_id(cell_ind) = 1;
            else
                sac_dir_id(cell_ind) = 0;
            end
            
            avg_STA_LFP = mean(STA_LFP);
            high_avg_STA_LFP = filtfilt(bhigh,ahigh,avg_STA_LFP);
            high_avg_STA_LFPd = mean(STA_LFPad);
            %high_avg_STA_LFPd = filtfilt(bhigh,ahigh,high_avg_STA_LFPd);
            
            rms = sqrt(mean(high_avg_STA_LFP.^2));
            [PKS,LOCS] = findpeaks(high_avg_STA_LFP,'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
            [trofs,trLOCS] = findpeaks(-high_avg_STA_LFP,'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
            
            if length(LOCS) == 1 %only want well detected ones
               trough_peak_aligned_STA(1,:) = trough_peak_aligned_STA(1,:)+avg_STA_LFP(LOCS-twin2:LOCS+twin2);
               tp_unit_count(1) = tp_unit_count(1)+1;
            end
            if length(trLOCS) == 1 %only want well detected ones
                trough_peak_aligned_STA(2,:) = trough_peak_aligned_STA(2,:)+avg_STA_LFP(trLOCS-twin2:trLOCS+twin2);
                tp_unit_count(2) = tp_unit_count(2)+1;
            end
            
            if spike_ind > 250
                avg_STA_LFP5 = NaN(5,twin*2+1);
                high_avg_STA_LFP5 = NaN(5,twin*2+1);
                quint_locs = NaN(1,5);
                quintile = floor((spike_ind-1)/5);
                for i = 1:5
                    avg_STA_LFP5(i,:) = mean(STA_LFP((i-1)*quintile+1:i*quintile,:));
                    high_avg_STA_LFP5(i,:) = filtfilt(bhigh,ahigh,avg_STA_LFP);
                    
                    if ~isempty(LOCS) && ~isempty(trLOCS) %% use the more prominent deviation for estimation
                        if PKS > -trofs
                            [pk,ql] = findpeaks(high_avg_STA_LFP5(i,:),'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
                            quint_locs(i) = ql(pk == max(pk));
                        else
                            [pk,ql] = findpeaks(-high_avg_STA_LFP5(i,:),'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
                            quint_locs(i) = ql(pk == max(pk));
                        end
                    elseif ~isempty(LOCS)
                        [pk,ql] = findpeaks(high_avg_STA_LFP5(i,:),'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
                        quint_locs(i) = ql(pk == max(pk));
                    elseif ~isempty(trLOCS)
                        [pk,ql] = findpeaks(-high_avg_STA_LFP5(i,:),'MinPeakHeight',3.7*rms,'MaxPeakWidth',15);
                        quint_locs(i) = ql(pk == max(pk));
                    end
                end
            end
            if ~all(isnan(quint_locs))
                peakvar =nanstd(quint_locs).^2;
            else
                peakvar = NaN;
            end
            
            %%
            
            recording_id(cell_ind) = monk*100+sess;
            recording_name{cell_ind} = [task_file(1:end-11) ' ' unit_stats{1,unit}];
            
            if ~isempty(trLOCS) && length(trLOCS) == 1 %%multiple peaks not good signal then
                min_max_lag(1,cell_ind) = trLOCS-twin-1;
            end
            if ~isempty(LOCS) &&  length(LOCS) == 1 %%multiple peaks not good signal then
                min_max_lag(2,cell_ind) = LOCS-twin-1;
            end
            
            within_sess_variance(cell_ind) = peakvar;
            
            cell_ind = cell_ind + 1;
            %%
%             figure
%             
%             subplot(2,3,1)
%             hold on
%             plot([-twin twin],[0 0],'k')
%             plot(-twin:twin,avg_STA_LFP)
%             if ~isempty(LOCS)
%                 plot(LOCS-twin-1,avg_STA_LFP(LOCS),'r*')
%             end
%             if ~isempty(trLOCS)
%                 plot(trLOCS-twin-1,avg_STA_LFP(trLOCS),'r*')
%             end
%             yl = ylim;
%             plot([0 0],[yl(1) yl(2)],'k--')
%             hold off
%             xlim([-twin twin])
%             xlabel('Time from Spike (ms)')
%             ylabel('LFP (uV)')
%             box off
%             title('Spike Triggered Average LFP (STA-LFP)')
%             
%             vals = STA_LFP(:);
%             vals(isnan(vals)) = [];
%             pct5 = prctile(vals,1);
%             pct95 = prctile(vals,99);
%             subplot(2,3,2)
%             imagesc(-twin:twin,1:size(STA_LFP,1),STA_LFP)
%             caxis([pct5 pct95])
%             xlabel('Time from Spike (ms)')
%             ylabel('Spike #')
%             title(['Channel #' num2str(unit_channel)])
%             box off
%             colormap('viridis')
%             title('Raw Spike Triggered LFP')
%             
%             subplot(2,3,5)
%             high_STA_LFP = STA_LFP;
%             for spk = 1:size(high_STA_LFP,1);
%                 high_STA_LFP(spk,:) = filtfilt(bhigh,ahigh,high_STA_LFP(spk,:));
%             end
%             imagesc(-twin:twin,1:size(high_STA_LFP,1),high_STA_LFP)
%             caxis([pct5 pct95])
%             xlabel('Time from Spike (ms)')
%             ylabel('Spike #')
%             title(['Channel #' num2str(unit_channel)])
%             box off
%             colormap('viridis')
%             title('High (20 Hz) Passed Spike Triggered LFP')
%             
%             subplot(2,3,4)
%             hold on
%             plot([-twin twin],[0 0],'k')
%             plot(-twin:twin,high_avg_STA_LFP)
%             if ~isempty(LOCS)
%                 plot(LOCS-twin-1,PKS,'r*')
%             end
%             if ~isempty(trLOCS)
%                 plot(trLOCS-twin-1,-trofs,'r*')
%             end
%             yl = ylim;
%             plot([0 0],[yl(1) yl(2)],'k--')
%             hold off
%             xlim([-twin twin])
%             xlabel('Time from Spike (ms)')
%             ylabel('LFP (uV)')
%             box off
%             title(['High (20 Hz) Passed STA-LFP: Max = ' num2str(LOCS-twin-1) ' ms, Min = ' num2str(trLOCS-twin-1) ' ms'])
%             
%             subplot(2,3,6)
%             hold on
%             plot([-twin twin],[0 0],'k')
%             plot(-twin:twin,avg_STA_LFP5')
%             yl = ylim;
%             plot([0 0],[yl(1) yl(2)],'k--')
%             hold off
%             xlim([-twin twin])
%             xlabel('Time from Spike (ms)')
%             ylabel('LFP (uV)')
%             box off
%             title(['STA-LFP by Qunitile: variance = ' num2str(peakvar) ' ms^2'])
%             
%             subplot(2,3,3)
%             hold on
%             plot([-twin twin],[0 0],'k')
%             ft = plot(-twin:twin,high_avg_STA_LFP);
%             %px = plot(-twin:twin,high_avg_STA_LFPd);
%             yl = ylim;
%             plot([0 0],[yl(1) yl(2)],'k--')
%             hold off
%             xlim([-twin twin])
%             xlabel('Time from Spike (ms)')
%             ylabel('LFP (uV)')
%             box off
%             %legend([ft,px],'Fieldtrip','Plexon')
%             title(['High (20 Hz) Passed STA-LFP: fieldtrip vs plexon import method'])
%             
%             subtitle([task_file(1:end-11) ' ' unit_stats{1,unit}])
%             %%
%             close
            %
        end
    end
end
%%
figure

min_lag_dist = min_max_lag(1,:);
min_lag_dist(isnan(min_lag_dist)) = [];
num_min = length(min_lag_dist);

subplot(2,2,1)
hist(min_lag_dist,25)
xlabel('Spike Time-Artifact (ms)')
ylabel('# of Units')
title(['Trough/Minimum Distribution, n = ' num2str(num_min), ', median = ' num2str(median(min_lag_dist)) ' ms'])
box off

max_lag_dist = min_max_lag(2,:);
max_lag_dist(isnan(max_lag_dist)) = [];
num_max = length(max_lag_dist);

subplot(2,2,2)
hist(max_lag_dist,25)
xlabel('Spike Time-Artifact (ms)')
ylabel('# of Units')
title(['Peak/Max Distribution, n = ' num2str(num_max) ', median = ' num2str(median(max_lag_dist)) ' ms'])
box off

win_sess_var = within_sess_variance;
win_sess_var(isnan(win_sess_var)) = [];
win_num = length(win_sess_var);

subplot(2,2,3)
hist(win_sess_var,25)
xlabel('Artifact-Artifact (ms^2)')
ylabel('# of Units')
title(['Within Unit Within Session Quantile Variance, n = ' num2str(win_num)])

rid = recording_id;
uniq_rid = unique(rid);
uniq_rid(isnan(uniq_rid)) = [];
rid_num_sess = length(uniq_rid);

within_sess_var = NaN(2,rid_num_sess);
for s = 1:length(uniq_rid)
    these_sess = find(rid == uniq_rid(s));
    these_lags1 =  min_max_lag(1,these_sess);
    within_sess_var(1,s) = nanstd(these_lags1).^2;
    
    these_lags2 =  min_max_lag(1,these_sess);
    within_sess_var(2,s) = nanstd(these_lags2).^2;
end

subplot(2,2,4)
hist(within_sess_var',25)
xlabel('Artifact-Artifact (ms^2)')
ylabel('# of Sessions')
title(['Across Units Within Session Variance, n = ' num2str(rid_num_sess)])
%%
average_record = laundry(average_image_ta_LFP);
figure
subplot(2,2,1)
imagesc(-twin:twin,1:size(average_image_ta_LFP,1),average_image_ta_LFP)
colormap('jet')
hold on
plot([0 0],[0 size(average_image_ta_LFP,1)],'w--')
hold off
xlabel('Time From Image Onset (ms)')
ylabel('LFP #')
title('Raw population data')
box off

[bhigh2,ahigh2] = butter(8,8/(Fs/2),'high');
[blow,alow] = butter(8,30/(Fs/2),'low');
filt_average_image_ta_LFP = NaN(size(average_image_ta_LFP));
trouch_time = NaN(1,size(average_image_ta_LFP,1));
for i = 1:size(filt_average_image_ta_LFP,1);
    filt_average_image_ta_LFP(i,:) = filtfilt(bhigh2,ahigh2,average_image_ta_LFP(i,:));
    strip = filt_average_image_ta_LFP(i,:);
    strip = filtfilt(blow,alow,strip);
    rms = sqrt(mean(strip.^2));
    [PKS,LOCS] = findpeaks(-strip,'MinPeakHeight',3*rms);
    if ~isempty(LOCS)
        trouch_time(i) = LOCS(max(PKS) == PKS);
    end
end

subplot(2,2,2)
imagesc(-twin:twin,1:size(filt_average_image_ta_LFP,1),filt_average_image_ta_LFP)
colormap('jet')
hold on
plot([0 0],[0 size(filt_average_image_ta_LFP,1)],'k--')
hold off
xlabel('Time From Image Onset (ms)')
ylabel('LFP #')
title('High (8 Hz) Passed population data')
box off


lag = min_max_lag(2,1:length(trouch_time));
nanlag = isnan(lag);
lag(nanlag) = [];
trt = trouch_time;
trt(nanlag) = [];
nantrt = isnan(trt);
trt(nantrt) = [];
lag(nantrt) = [];
[r,p] = corrcoef(trt,lag);
P = polyfit(lag,trt,1);


subplot(2,2,3)
lag = min_max_lag(2,1:length(trouch_time));
plot(trouch_time-twin,lag,'.k')
xlabel('Trouch Relative to Image Onset (ms)')
ylabel('Spike-LFP artifact peak (ms)')
box off
title(['Corr btwn Image Onset Lag and Spike Artifact Lag: m = ' num2str(P(1),3) ', r  = ' num2str(r(2),3) ', p = ' num2str(p(2),3)])
axis equal

sorted_filt_average_image_ta_LFP = filt_average_image_ta_LFP;
[s,si] = sort(min_max_lag(2,1:length(trouch_time)));
sorted_filt_average_image_ta_LFP = sorted_filt_average_image_ta_LFP(si,:);
nanind = find(isnan(s));
nanind = nanind(1);

subplot(2,2,4)
imagesc(-twin:twin,1:size(sorted_filt_average_image_ta_LFP,1),sorted_filt_average_image_ta_LFP)
hold on
plot([0 0],[0 size(sorted_filt_average_image_ta_LFP,1)],'k--')
plot([-twin twin],[nanind nanind],'k--')
hold off
colormap('jet')
hold on
plot([0 0],[0 size(sorted_filt_average_image_ta_LFP,1)],'w--')
hold off
xlabel('Time From Image Onset (ms)')
ylabel('LFP #')
title('Sorted by Spike-LFP lag & High (8 Hz) Passed population data')
box off
%%
figure
subplot(1,2,1)
plot(-twin2:twin2,trough_peak_aligned_STA(1,:)/tp_unit_count(1));
hold on
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
box off
xlabel('Time from Peak (ms)')
ylabel('LFP (uV)')
title(['Peak Aligned, n = ' num2str(tp_unit_count(1))])

subplot(1,2,2)
plot(-twin2:twin2,trough_peak_aligned_STA(2,:)/tp_unit_count(2));
hold on
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
box off
xlabel('Time from Trough (ms)')
ylabel('LFP (uV)')
title(['Trough Aligned, n = ' num2str(tp_unit_count(2))])

subtitle('Average Spike-LFP Artifact')


%% Need to generate plot
start_vals = (first_last_spike(1,:)-analog_start_end_ts(1,:));
start_vals(start_vals > 100) = 100;

end_vals = (first_last_spike(2,:)-analog_start_end_ts(2,:));
end_vals(end_vals <  -100) = -100;

figure
subplot(1,2,1)
plot(min_max_lag(2,:),start_vals,'.k')
xlabel('Spike-LFP Artifact Lag (ms)')
ylabel('1st Thershold Crossing-1st Analog Time Stamp (ms)')
box off
axis square
title('Spike Channel Start Time-Analog Channel Start time')

subplot(1,2,2)
plot(min_max_lag(2,:),end_vals,'.k')
xlabel('Spike-LFP Artifact Lag (ms)')
ylabel('Last Thershold Crossing-Last Analog Time Stamp (ms)')
box off
axis square
title('Spike Channel End Time-Analog Channel End time')