function make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data,task)
% written by Seth Konig August 2014. Updated January 2016 by SDK
% 1) Plots spike trains aligned to trial start. Below is the PSTH.
% 2) Plots waveform divided into 1/4 of the task at time
% 3) Plots firing rate over time/trials
% 4) defines which trials should be used for data analysis if firing rate
% is not stable

%aligns data to item 1 on for sequence task
%aligns data to image onset for image task
%aligns data to dot onset for cvtnew task
%also aligns data to ITI since there seems to be a lot of reward/ITI neurons

% don't want to waste processing power on units that are a) poorly isolated,
% b) don't have enough spikes, and c) aren't stable for long enough to do
% any analysis with sufficient confidence. But should plot the rasters
% nevertheless.
%!!!WARNING CODE WILL SET VIABLE TRILAS TO NULL FOR THESE UNITS TO SAVE PROCESSING TIME!!!!
%
% Code rechecked for bugs August, 2016

figure_dir = [figure_dir 'Rasters\'];
too_few_folder = [figure_dir '\Rasters Too Few Trials or Spikes\']; %where to put rasters see line 47
PoorIsolationFolder = [figure_dir '\PoorIsolation\']; %where to put poorly isolated units
MultiUnit_folder = [figure_dir '\MultiUnit\'];


reward_code = 3;
ITI_code = 15;
dot_on_code = 25;%for cvtnew
trial_start_code = 23; %image on and item 1 on
binsize = 35;%ms probably want 25-100 ms
trialbinsize = 5;%averging data by trial number,
%used 6 in the past 5 makes more sennse since 20 famiarization trials and
%firing rates can vary drastically across tasks

%for neurons with too few spikes or too few stable trials. Hard to properly
%isolate but also over isolation period has a very low average firing rate.
% Prior experinece suggest these neurons typically don't fire except
% sporadically in burst or on specifc fixations. May be useful to analyze
% in the future.
poor_isolation_threshold = 2;%0,1,2 on cluster cutting quality
confidence_threshold = 0.79;% percent condifence that I'm willing to take anything over 80%

%get important task related data
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~,~,comments]=...
    get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end

%these are the absolute minimum data required to do data analysis may want
%to be more strigent later but not worth doing analysis (especially
%shuffling) on such few trials for these neurons
if strcmpi(task,'cvtnew')
    minimum_trials_1 = 100;
    load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','data','hdr','waveforms');
    minimum_spikes = 75; %slightly less than ListSQ because shorter period but often more stable too
else %for listSQ
    if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
        minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
        minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
    else
        minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
        minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
    end
    minimum_spikes = 100; %currently setup for all trials excludes before and after task was running
    load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','data','hdr','fixationstats','waveforms');
end

%get unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality,comments);

if num_units == 0;
    return;
end

%---Import Channel Parameters---%
%get channel parameters like gain and threshold. Max voltage is 3000 uV
%have to import using .plx file not .nex file it's pretty quick anyway
[~,names] = plx_chan_names([data_dir task_file(1:end-3) 'plx']);
chan_num = NaN(1,4);
for i = 1:4
    chan_num(i) = str2double(names(i,6));
end
clear names
[~,gains] = plx_chan_gains([data_dir task_file(1:end-3) 'plx']);
gains = gains(1:4);
[~,thresholds] = plx_chan_thresholds([data_dir task_file(1:end-3) 'plx']);
thresholds = thresholds(1:4);

num_trials = length(cfg.trl);
switch task
    case 'ListSQ'
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        
        valid_trials = NaN(2,num_units);%trials we want to keep based on firing rate
        stability_attribute = ones(1,num_units); %1 stable and using, 2 insufficient number of spikes,
        %3 not stable for sufficient trials, 4 poorly isolated/not confident unit is real
        for unit = 1:num_units
            
            type = NaN(1,length(cfg.trl));
            allspikes = NaN(length(cfg.trl),7000);
            allspikesITI = NaN(length(cfg.trl),1750);
            allsaccadespikes = NaN(length(cfg.trl),1000); %looking at eye related activity mostly for low firing rate neurons
            
            gain = gains(str2double(unit_stats{1,unit}(6)));
            waveforms{1,unit} = gain*waveforms{1,unit}; %conver to uV
            avg_waveform_amplitude = max(waveforms{1,unit})-min(waveforms{1,unit});
            threshold = thresholds(str2double(unit_stats{1,unit}(6)));
            
            %calculate average waveform in terms of 10ths of session that
            %you have the neuron for
            avg_waveform = cell(1,10);
            segments = floor(size(waveforms{1,unit},2)/10);
            for i = 1:10;
                avg_waveform{i}= mean(waveforms{1,unit}(:,(segments*(i-1)+1):(segments*i))');
            end
            
            
            spike_count = length(avg_waveform_amplitude);
            trial_averaged_amplitude = NaN(1,num_trials);
            trial_num = NaN(1,num_trials);
            for t = 1:num_trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITI_code); %ms sample index
                event = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code)-trial_start; %event start within trial time
                
                if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) % for the sequence trials
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 1;
                    else
                        continue;
                    end
                elseif  itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2)
                    if  any(cfg.trl(t).allval == reward_code)
                        type(t) = 2;
                    else
                        continue;
                    end
                else
                    if  any(cfg.trl(t).allval == trial_start_code) %if image was shown
                        type(t) = 3;
                    else
                        continue;
                    end
                end
                
                trial_num(t) = t;
                
                %locked to the ITI
                spikes =data(unit).values{t};
                allspikesITI(t,:) = spikes(1:1750);
                
                %locked to main event                
                waveform_ind = find(spikes)+trial_start-1;
                ind = NaN(1,length(waveform_ind));
                for ii = 1:length(waveform_ind)
                    %give 1 ms quantization error buffer 1 ms for spike
                    %times + 1 ms for trial start = 2 ms
                    temp = find(waveforms{2,unit} <= waveform_ind(ii)+2 & waveforms{2,unit} >= waveform_ind(ii)-2);
                    if ~isempty(temp) %not sure why this happens but some trials are padded by NaNs so timing could be off within certain trials
                        temp = temp(1); %only take the first one in case several spikes occur within a short period of time
                        ind(ii) = temp;
                    end
                end
                if ~isempty(ind)
                    ind(isnan(ind)) = [];
                    trial_averaged_amplitude(t) = nanmean(avg_waveform_amplitude(ind)); %grab waveforms from trial and average
                else
                    trial_averaged_amplitude(t) = NaN;
                end
                
                spikes = find(spikes);
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 7000) = [];
                spks = zeros(1,7000);
                spks(spikes) = 1;
                allspikes(t,:) = spks;
                
                
                if type(t) == 3 %if image type look at saccade aligned activity
                    saccadetimes = fixationstats{t}.saccadetimes;
                    invalid= find(saccadetimes(1,:) < event+500); %ignore first 500 ms
                    saccadetimes(:,invalid) = [];
                    
                    saccade_algined = zeros(1,1000);
                    
                    %going to average across all saccades mostly for low firing rate neurons
                    for s = 1:size(saccadetimes,2)
                        spks = spikes-saccadetimes(1,s);
                        spks(spks < -499) =[];
                        spks(spks > 500) = [];
                        saccade_algined(spks+500) = 1; %not an error only care when spikes not amount
                    end
                    allsaccadespikes(t,:) = saccade_algined;
                end
            end
            
            %don't want to clean these up yet since trial # is important
            %determine spike per groups of trials to determine if firing rate is approximately stable over time
            spikespertrial = bin1(allspikes',trialbinsize)'./trialbinsize;
            ITIspikespertrial = bin1(allspikesITI',trialbinsize)'./trialbinsize;
            allspikespertrial = spikespertrial+ITIspikespertrial; %not mutually exclusive but doens't matter for this
            trial_averaged_amplitude = bin1(trial_averaged_amplitude,trialbinsize,'mean');
            %have same number as bins as trial data though time of spikes
            %may be slightly off
            
            %cleanup by remove unsuccessful trial #'s
            type = laundry(type);
            allspikes = laundry(allspikes,1);
            allspikesITI = laundry(allspikesITI,1);
            allsaccadespikes = laundry(allsaccadespikes,1);
            trial_averaged_amplitude(trial_averaged_amplitude == 0) = NaN;
            %replace 0's with NaNs didn't want to erase the data before
            trial_num = laundry(trial_num);
            
            title_str = unit_stats{1,unit};
            if  multiunit(unit)
                title_str = ['Multiunit '  title_str];
            end
            title_str = [title_str ' Unit Confidence ' num2str(100*unit_stats{2,unit}) '% ' ...
                'Cutting Quality ' num2str(unit_stats{3,unit})];
            
            if spike_count < minimum_spikes %too few spikes
                listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
                    trial_averaged_amplitude,trial_num,avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit})
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 2;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            elseif unit_stats{3,unit} <= poor_isolation_threshold || unit_stats{2,unit} <= confidence_threshold
                %not confident it's a unit or good isolation
                listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
                    trial_averaged_amplitude,trial_num,avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit})
                save_and_close_fig(PoorIsolationFolder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 4;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            end
            
            listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
                ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
                trial_averaged_amplitude,trial_num,avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit})
            
            %do manual check of which trials to take. Humans are just
            %really good at seeing pattens....
            try
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            catch
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            end
            if isnumeric(reply)
                [start_end] = listsqTrialCutoffplots(reply,allspikes,type,allspikesITI,allspikespertrial,...
                    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
                    trial_averaged_amplitude,trial_num,avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit});
            else
                start_end = [NaN NaN]; %then take all trials
            end
            
            %determine if the number of trials is actually sufficient to
            %use for data analysis
            if start_end(2) > length(cfg.trl) %human error
               start_end(2) = length(cfg.trl); %fix the stupid humans error 
            end
            if any(~isnan(start_end)) && all(start_end ~= 0); %if all NaNs take all trials
                if isnan(start_end(2))
                    max_trial = cfg.trl(t).cnd-1000; %get the last condition number
                    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
                    if max_trial >  minimum_trials_1 %if last trial was greater than the end of the first novel/repeat block
                        if max_trial-min_trial < minimum_trials_2 %not enough trials didn't get through 1 novel/repeat block
                            start_end = [0 0];
                        end
                    else %didn't even get through first novel/repeat block
                        start_end = [0 0];
                    end
                elseif isnan(start_end(1)) %then take first trial to max_trial
                    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
                    if max_trial < minimum_trials_1 %didn't finish fam block + at least 1 novel/repeat block
                        start_end = [0 0];
                    end
                else %both trials numeric
                    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the last condition number
                    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
                    if max_trial < minimum_trials_1 %didn't finish fam block + at least 1 novel/repeat block
                        start_end = [0 0];
                    elseif max_trial-min_trial < minimum_trials_2 %not enough trials didn't get through 1 novel/repeat block
                        start_end = [0 0];
                    end
                end
            end
            
            valid_trials(:,unit) = start_end';
            
            subtitle(title_str);
            
            if all(valid_trials(:,unit) == 0)%not enough trials to do data analysis on
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                stability_attribute(unit) = 3;
            else
                if multiunit(unit) == 1 %then treat as multiunit
                    save_and_close_fig(MultiUnit_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                else
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                end
            end
        end
    case {'cvtnew','CVTNEW'}
        
        valid_trials = NaN(2,num_units);%trials we want to keep based on firing rate
        stability_attribute = ones(1,num_units); %1 stable and using, 2 insufficient number of spikes,
        %3 not stable for sufficient trials, 4 poorly isolated/not confident unit is real
        for unit = 1:num_units
            allspikes = NaN(length(cfg.trl),3500);
            allspikesITI = NaN(length(cfg.trl),1500);
            rewardspikes = NaN(length(cfg.trl),1500);
            
            gain = gains(str2double(unit_stats{1,unit}(6)));
            waveforms{1,unit} = gain*waveforms{1,unit}; %conver to uV
            avg_waveform_amplitude = max(waveforms{1,unit})-min(waveforms{1,unit});
            threshold = thresholds(str2double(unit_stats{1,unit}(6)));
            
            %calculate average waveform in terms of 10ths of session that
            %you have the neuron for
            avg_waveform = cell(1,10);
            segments = floor(size(waveforms{1,unit},2)/10);
            for i = 1:10;
                avg_waveform{i}= mean(waveforms{1,unit}(:,(segments*(i-1)+1):(segments*i))');
            end
            
            spike_count = length(avg_waveform_amplitude);
            trial_num = NaN(1,num_trials);
            for t = 1:num_trials
                if ~any(cfg.trl(t).allval == reward_code) %successful trial
                    continue
                end
                trial_num(t) = t;
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITI_code);
                event = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start;
                reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start;
                reward = reward(1);
                
                spikes =data(unit).values{t};
                
                
                waveform_ind = find(spikes)+trial_start-1;
                ind = NaN(1,length(waveform_ind));
                for ii = 1:length(waveform_ind)
                    %give 1 ms quantization error buffer 1 ms for spike
                    %times + 1 ms for trial start = 2 ms
                    temp = find(waveforms{2,unit} <= waveform_ind(ii)+2 & waveforms{2,unit} >= waveform_ind(ii)-2);
                    temp = temp(1); %only take the first one in case several spikes occur within a short period of time
                    ind(ii) = temp;
                end
                if ~isempty(ind)
                    trial_averaged_amplitude(t) = nanmean(avg_waveform_amplitude(ind)); %grab waveforms from trial and average
                else
                    trial_averaged_amplitude(t) = NaN;
                end
                
                %locked to the ITI
                allspikesITI(t,:) = spikes(1:1500);
                
                %locked to reward
                rewardspikes(t,:) = spikes(reward-500:reward+999);
                
                %locked to main event
                spikes = find(spikes);
                spikes = spikes-event;
                spikes(spikes < 1) = [];
                spikes(spikes > 3500) = [];
                temp = zeros(1,3500);
                temp(spikes) = 1;
                allspikes(t,:) = temp;
            end
            
            
            %determine spike per groups of trials to determine if firing rate is approximately stable over time
            spikespertrial = bin1(allspikes',trialbinsize)'./trialbinsize;
            ITIspikespertrial = bin1(allspikesITI',trialbinsize)'./trialbinsize;
            allspikespertrial = spikespertrial+ITIspikespertrial; %not mutually exclusive but doens't matter for this
            trial_averaged_amplitude = bin1(trial_averaged_amplitude,trialbinsize,'mean');
            
            %cleanup by remove unsuccessful trial #'s
            allspikes = laundry(allspikes,1);
            allspikesITI = laundry(allspikesITI,1);
            trial_averaged_amplitude(trial_averaged_amplitude == 0) = NaN;
            rewardspikes = laundry(rewardspikes,1);
            %replace 0's with NaNs didn't want to erase the data before
            trial_num = laundry(trial_num);
            
            title_str = unit_stats{1,unit};
            if  multiunit(unit)
                title_str = ['Multiunit '  title_str];
            end
            title_str = [title_str ' Unit Confidence ' num2str(100*unit_stats{2,unit}) '% ' ...
                'Cutting Quality ' num2str(unit_stats{3,unit})];
            
            if spike_count < minimum_spikes %too few spikes
                
                cvtnewplot(allspikes,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,rewardspikes,trialbinsize,binsize,trial_averaged_amplitude,...
                    avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit},trial_num);
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 2;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            elseif unit_stats{3,unit} <= poor_isolation_threshold || unit_stats{2,unit} <= confidence_threshold
                %not confident it's a unit or good isolation
                cvtnewplot(allspikes,allspikesITI,allspikespertrial,spikespertrial,...
                    ITIspikespertrial,rewardspikes,trialbinsize,binsize,trial_averaged_amplitude,...
                    avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit},trial_num);
                save_and_close_fig(PoorIsolationFolder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                
                stability_attribute(unit) = 4;
                valid_trials(:,unit) = 0; %don't want to process further so don't analyze any trials for this unit
                continue %make raster then go on to next neuron
            end
            
            
            cvtnewplot(allspikes,allspikesITI,allspikespertrial,spikespertrial,...
                ITIspikespertrial,rewardspikes,trialbinsize,binsize,trial_averaged_amplitude,...
                avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit},trial_num);
            
            
            %do manual check of which trials to take. Humans are just
            %really good at seeing pattens....
            try
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            catch
                reply = input(['You can keep all trials (max trial #' num2str(num_trials) '). Is this Ok? [Y/#s]: \n' ...
                    'If NO please state [trial start and trial end].']);
            end
            if isnumeric(reply)
                [start_end] = cvtnewTrialCutoffplots(reply,allspikes,allspikesITI,rewardspikes,allspikespertrial,...
                    spikespertrial,ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,...
                    avg_waveform,waveforms(:,unit),threshold,title_str,unit_stats{4,unit},trial_num);
            else
                start_end = [NaN NaN]; %then take all trials
            end
            valid_trials(:,unit) = start_end';
            
            
            %determine if the number of trials is actually sufficient to
            %use for data analysis
            if any(~isnan(start_end)) && all(start_end ~= 0); %if all NaNs take all trials
                if isnan(start_end(2))
                    max_trial = cfg.trl(t).cnd-1000; %get the last condition number
                else
                    max_trial = trial_num(end);
                end
                if isnan(start_end(1)) %then take first trial to max_trial
                    min_trial = 1;
                else
                    min_trial = valid_trials(1,unit);
                end
                
                if max_trial-min_trial < minimum_trials_1
                    valid_trials(:,unit) = 0;
                end
            end
            
            
            subtitle(title_str);
            if all(valid_trials(:,unit) == 0)%not enough trials to do data analysis on
                save_and_close_fig(too_few_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                stability_attribute(unit) = 3;
            else
                if multiunit(unit) == 1 %then treat as multiunit
                    save_and_close_fig(MultiUnit_folder,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                else
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_raster']);
                end
            end
            
        end
end

%add the valid trials variable to preprocess file
save([data_dir task_file(1:end-11) '-preprocessed.mat'],'-append','valid_trials','stability_attribute')
end

function listsqplot(allspikes,type,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
    trial_averaged_amplitude,trial_num,avg_waveform,waveforms,threshold,title_str,comments)
%plot the listsq

if ~isempty(findall(0,'Type','Figure'))
    g = gcf;
    if g.Number == 101;
        close
    end
end

% Rasters from Main Event
maxtime = 5000*ones(1,3);

figure(101);

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
pause(0.5) %give time for plot to reach final size

%Raster from ITI start
subplot(3,3,1)
[trial,time] = find(allspikesITI == 1);
if ~isempty(trial)
    plot(time,trial_num(trial),'.k')
    ylim([0 max(trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from ITI start (ms)')

%Raster from Sequence Trial Start
sub_trial_num =	trial_num(find(type == 1 | type == 2 )); %sequence items
subplot(3,3,2)
[trial,time] = find(allspikes(type == 1 | type == 2 ,:)== 1);
plot(time,sub_trial_num(trial),'.k')
if ~isempty(trial)
    ylim([0 max(sub_trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from Sequence Trial start (ms)')
if ~isempty(time)
    if max(time) > 5000
        xlim([0 max(time)])
    else
        xlim([0 max(time)])
    end
end

%Time from Image Trial Start
sub_trial_num =	trial_num(find(type == 3)); %sequence items
subplot(3,3,3)
[trial,time] = find(allspikes(type == 3,:)== 1);
plot(time,sub_trial_num(trial),'.k')
if ~isempty(trial)
    ylim([0 max(sub_trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from Image Onset (ms)')
if ~isempty(time)
    if max(time) > 5000
        xlim([0 max(time)])
    else
        xlim([0 max(time)])
    end
end

%just raster from saccade aligned activtiy. Any spike aligned to a saccade
%mostly for low firing rate neurons to determine stability for these when
%not obvious
subplot(3,3,4)
[trial,time] = find(allsaccadespikes == 1);
if ~isempty(trial)
    plot(time-500,sub_trial_num(trial),'.k')
    ylim([0 max(sub_trial_num(trial))+5])
end
ylabel('Trial #')
xlabel('Time from Saccade Start (ms)')

subplot(3,3,5)
hold all
for i = 1:length(avg_waveform)
    plot(avg_waveform{i})
end
plot([1 32],[threshold threshold],'k--')
hold off
legend(cellstr(num2str((1:10)')))
ylabel('Avg Waveform (uV)')
xlim([1 32])
ylim([-3000 3000]) %max voltage is 3000 uV

subplot(3,3,6)
plot(max(waveforms{1})-min(waveforms{1}),waveforms{2},'k.')
xlabel('Amplitude (uV)')
ylabel('"Time Stamp"')
axis tight


allspikespertrial(allspikespertrial == 0) = NaN;
subplot(3,3,7:8)
hold on
b = bar([spikespertrial ITIspikespertrial],'stacked');
set(b,'edgecolor','none','FaceAlpha',0.5)
if ~isempty(trial)
    xlim([0 find(~isnan(allspikespertrial),1,'last')+1]);
end
xl = xlim;
plot([xl(1) xl(2)],[nanmedian(allspikespertrial) nanmedian(allspikespertrial)],'k-','linewidth',5)
plot([xl(1) xl(2)],[nanmedian(allspikespertrial)-nanstd(allspikespertrial) nanmedian(allspikespertrial)-nanstd(allspikespertrial)],'k--')
plot([xl(1) xl(2)],[nanmedian(allspikespertrial)+nanstd(allspikespertrial) nanmedian(allspikespertrial)+nanstd(allspikespertrial)],'k--')
yl = ylim;
ylim([0 yl(2)]);

%scale waveform amplitude so sits nicely on plot
trial_averaged_amplitude = trial_averaged_amplitude-nanmean(trial_averaged_amplitude); %zero
trial_averaged_amplitude = trial_averaged_amplitude/(max(abs(trial_averaged_amplitude))); %scale to 1
trial_averaged_amplitude = yl(2)/3*trial_averaged_amplitude; %rescale scale
avg_waveform_ampltiude = trial_averaged_amplitude+yl(2)/2;%set level to median
plot(avg_waveform_ampltiude,'r')

hold off
xlabel(['Groups of Trials (' num2str(trialbinsize) 'trials/group)'])
ylabel('Average Spikes per trial')

subplot(3,3,9)
len = length(comments);
rows = floor(len/40);
if rows > 1
    cmnts = [];
    for r = 1:rows
        rind = 40*(r-1)+1:40*r;
        if r == rows
        cmnts = [cmnts '\n' comments(rind(1):end)];
        else
            cmnts = [cmnts '\n' comments(rind)];
        end
    end
    text(0.5,0.5,sprintf(cmnts))
else
    if ~isnan(comments) %if NaN no comments
        text(0.5,0.5,comments)
    end
end
title('Comments')
axis off


subtitle(title_str);

end

function  cvtnewplot(allspikes,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,rewardspikes,trialbinsize,binsize,trial_averaged_amplitude,...
    avg_waveform,waveforms,threshold,title_str,comments,trial_num)
%plot the cvtnew

if ~isempty(findall(0,'Type','Figure'))
    g = gcf;
    if g.Number == 101;
        close
    end
end


figure(101);

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)]);
pause(0.5) %give time for plot to reach final size


% Rasters from Main Event
subplot(3,3,1)
[trial,time] = find(allspikesITI == 1);
if ~isempty(trial)
    plot(time,trial_num(trial),'.k')
    ylim([0 max(trial_num(trial))+5])
end
ylim([0 max(trial)+5])
ylabel('Trial #')
xlim([0 1500])
xlabel('Time from ITI start (ms)')
box off

subplot(3,3,2)
[trial,time] = find(allspikes == 1);
if ~isempty(trial)
    plot(time,trial_num(trial),'.k')
    ylim([0 max(trial_num(trial))+5])
end
ylabel('Trial #')
xlim([0 3500])
xlabel('Time from Dot on (ms)')
box off

subplot(3,3,3)
[trial,time] = find(rewardspikes == 1);
if ~isempty(trial)
    plot(time,trial_num(trial),'.k')
    ylim([0 max(trial_num(trial))+5])
end
ylabel('Trial #')
xlim([0 1500])
xlabel('Time from Reward start (ms)')
box off

subplot(3,3,[7 8])
hold on
b = bar([spikespertrial ITIspikespertrial],'stacked');
if ~isempty(trial)
    xlim([0 find(~isnan(spikespertrial),1,'last')+1]);
end
set(b,'edgecolor','none','FaceAlpha',0.5)
xl = xlim;
plot([xl(1) xl(2)],[median(allspikespertrial) median(allspikespertrial)],'k-','linewidth',5)
plot([xl(1) xl(2)],[median(allspikespertrial)-std(allspikespertrial) median(allspikespertrial)-std(allspikespertrial)],'k--')
plot([xl(1) xl(2)],[median(allspikespertrial)+std(allspikespertrial) median(allspikespertrial)+std(allspikespertrial)],'k--')
yl = ylim;
ylim([0 yl(2)]);

%scale waveform amplitude so sits nicely on plot
trial_averaged_amplitude = trial_averaged_amplitude-mean(trial_averaged_amplitude); %zero
trial_averaged_amplitude = yl(2)/2*trial_averaged_amplitude/max(abs(trial_averaged_amplitude)); %scale
avg_waveform_ampltiude = trial_averaged_amplitude+yl(2)/2;%set level to median
plot(avg_waveform_ampltiude,'r','linewidth',3)

hold off
xlabel(['Groups of Trials (' num2str(trialbinsize) 'trials/group)'])
ylabel('Average Spikes per trial')


subplot(3,3,5)
hold all
for i = 1:length(avg_waveform)
    plot(avg_waveform{i})
end
plot([1 32],[threshold threshold],'k--')
hold off
legend(cellstr(num2str((1:10)')))
ylabel('Avg Waveform (uV)')
xlim([1 32])
ylim([-3000 3000]) %max voltage is 3000 uV

subplot(3,3,6)
plot(max(waveforms{1})-min(waveforms{1}),waveforms{2},'k.')
xlabel('Amplitude (uV)')
ylabel('"Time Stamp"')
axis tight
box off


subplot(3,3,9)
len = length(comments);
rows = floor(len/40);
if rows > 1
    cmnts = [];
    for r = 1:rows
        rind = 40*(r-1)+1:40*r;
        if r == rows
        cmnts = [cmnts '\n' comments(rind(1):end)];
        else
            cmnts = [cmnts '\n' comments(rind)];
        end
    end
    text(0.5,0.5,sprintf(cmnts))
else
    if ~isnan(comments) %if NaN no comments
        text(0.5,0.5,comments)
    end
end
title('Comments')
axis off


subtitle(title_str);

end

function [start_end] = listsqTrialCutoffplots(reply,allspikes,type,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
    trial_averaged_amplitude,trial_num,avg_waveform,waveforms,threshold,title_str,comments)
%plot the cutoffs from replys
while isnumeric(reply)
    listsqplot(allspikes,type,allspikesITI,allspikespertrial,spikespertrial,...
        ITIspikespertrial,allsaccadespikes,trialbinsize,binsize,...
        trial_averaged_amplitude,trial_num,avg_waveform,waveforms,threshold,title_str,comments)
    
    start_end = reply;
    
    if length(reply) ~= 2
       reply = [0 0]; 
    end
    
    %for plotting and visualization should put a line down
    nano = 0;
    if isnan(reply(1));
        reply(1) = 0;
    end
    if isnan(reply(2))
        nano = 1;
    end
    for sb = [1:4 7];
        if sb == 7
            subplot(3,3,7:8)
            yl = ylim;
            hold on
            plot([reply(1)/trialbinsize reply(1)/trialbinsize],[0 yl(2)],'r');
            plot([reply(2)/trialbinsize reply(2)/trialbinsize],[0 yl(2)],'r');
            hold off
        else
            subplot(3,3,sb)
            xl = xlim;
            if nano
                yl = ylim;
                hold on
                plot([0 xl(2)],[reply(1) reply(1)],'r');
                plot([0 xl(2)],[yl(2) yl(2)],'r');
                hold off
            else
                hold on
                plot([0 xl(2)],[reply(1) reply(1)],'r');
                plot([0 xl(2)],[reply(2) reply(2)],'r');
                hold off
            end
        end
    end
    subtitle(title_str);
    
    try
        reply = input(['You can keep these trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    catch
        reply = input(['You can keep these trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    end

    
end
end

function [start_end] = cvtnewTrialCutoffplots(reply,allspikes,allspikesITI,rewardspikes,allspikespertrial,...
    spikespertrial,ITIspikespertrial,trialbinsize,binsize,trial_averaged_amplitude,avg_waveform,...
    waveforms,threshold,title_str,comments,trial_num)
%plot the cutoffs from replys
while isnumeric(reply)
    
    cvtnewplot(allspikes,allspikesITI,allspikespertrial,...
    spikespertrial,ITIspikespertrial,rewardspikes,trialbinsize,binsize,trial_averaged_amplitude,...
    avg_waveform,waveforms,threshold,title_str,comments,trial_num)

    start_end = reply;        
    if length(reply) ~= 2
       reply = [0 0]; 
    end
    
    
    %for plotting and visualization should put a line down
    nano = 0;
    if isnan(reply(1));
        reply(1) = 0;
    end
    if isnan(reply(2))
        nano = 1;
    end
    
    subplot(3,3,[7 8])
    yl = ylim;
    hold on
    if nano
        plot([reply(1)/trialbinsize reply(1)/trialbinsize],[0 yl(2)],'r');
        plot([size(allspikespertrial,1) size(allspikespertrial,1)],[0 yl(2)],'r');
    else
        plot([reply(1)/trialbinsize reply(1)/trialbinsize],[0 yl(2)],'r');
        plot([reply(2)/trialbinsize reply(2)/trialbinsize],[0 yl(2)],'r');
    end
    hold off
    
    subplot(3,3,1)
    xl = xlim;
    if nano
        yl = ylim;
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[yl(2) yl(2)],'r');
        hold off
    else
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[reply(2) reply(2)],'r');
        hold off
    end
    
    subplot(3,3,2)
    xl = xlim;
    if nano
        yl = ylim;
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[yl(2) yl(2)],'r');
        hold off
    else
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[reply(2) reply(2)],'r');
        hold off
    end
    
    subplot(3,3,3)
    xl = xlim;
    if nano
        yl = ylim;
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[yl(2) yl(2)],'r');
        hold off
    else
        hold on
        plot([0 xl(2)],[reply(1) reply(1)],'r');
        plot([0 xl(2)],[reply(2) reply(2)],'r');
        hold off
    end
    
    try
        reply = input(['You can keep all trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    catch
        reply = input(['You can keep all trials. Is this Ok? [Y/#s]: \n' ...
            'If NO please state [trial start and trial end].']);
    end
end

end

