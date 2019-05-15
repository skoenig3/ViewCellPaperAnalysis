function  Place_Cell_FreqModulation(data_dir,figure_dir,session_data,task)
%written by Seth Konig 7/23/6
%code only analyzes units with significant (95+ percentile) skaggs information
%score and Miriam's spatial stability score. Code combines analysis across
%various view cells codes to determine if neruons are theta modulated or
%possibly modulated by other frequncies. Code looks at all spikes, spikes
%on during fixations inside the view field vs outside the view field, and
%spikes following fixation in the view field vs preceding the fixation in
%the view field. Code analyzes frequency modulated in 3 ways: 1) STA LFP,
%2) spike train autocorrelatin, and 3) spike field coherence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = 500;% here how much after image onset to ignore data
view_cell_dir = [figure_dir 'View Cells\FreqModulation\'];
imageX = 800;
imageY = 600;
img_on_code = 23;
img_off_code = 24;
ITIstart_code = 15;
min_bin_dur = 0.100; %minimum of 100 ms in each bin to use so no outlier


smval = [60 10]; %1 is for fixation locked firing rate, 2 is for autocorr
minimum_fix_duration = 100;%miniumum fixation duration to look at data for
smval = [60 10];%gaussian 1/2 width for smoothing

%spike-feld coherence variables
buffer = 2048;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;

task_file = get_task_data(session_data,task);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
load([data_dir task_file(1:end-11) '-preprocessed.mat']);
absolute_fixationstats = fixationstats;
absolute_cfg = cfg;

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~,lfp_quality]...
    = get_task_data(session_data,task);
[~,~,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);

%only need right now when formatting is a little weird while resorting
unit_names = cell(1,num_units);
for unit = 1:num_units
    unit_names{unit}= hdr.label{data(unit).whichchannel};
end

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
    (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%NaNs are for start and end trials otherwise cut
num_trials = length(cfg.trl);
valid_trials(1,isnan(valid_trials(1,:))) = 1;
valid_trials(2,isnan(valid_trials(2,:))) = num_trials;

%these are the absolute minimum data required to do data analysis may want
%to be more strigent later but not worth doing analysis (especially
%shuffling) on such few trials for these neurons
if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
    minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
else
    minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
    minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
end


for unit = 1:size(valid_trials,2)
    %---determine number of blocks unit was stable for
    start_end = valid_trials(:,unit);
    if isnan(start_end(1))
        start_end(1) = 1;
    end
    if isnan(start_end(2))
        start_end(2) = length(cfg.trl);
    end
    start_end(start_end == 0) = 1;
    min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
    max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
    
    if min_trial < 22 %includes fam block
        min_trial = 22; %remove then count from there
    end
    num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
    
    if num_blks < 2
        valid_trials(:,unit) = 0;
    end
end


%filter parameters
filter_size = filter_width*10;
H = fspecial('gaussian',filter_size,filter_width);

%---Pre-allocate space for Fixations Inside vs Outside Field Analysis---%
list_fixation_locked_firing = cell(1,num_units); %firing rate locked to fixations
in_out = cell(1,num_units); %firing rate for fixations inside (1) or outside (0) of view filed
area = NaN(1,num_units); %area of place field
list_STA_LFP = cell(1,num_units); %Spike Triggered Average (STA) LFP
spike_autocorr = cell(2,num_units);%row 1 inside field, row 2 outside field ignores fixations without spikes
all_low_statSts = cell(1,num_units); %low frequency spike_field_coherence 
all_high_statSts = cell(1,num_units); %high frequency spike_field_coherence 

for unit = 1:num_units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Unit Passes 95% for both Skaggs and Stability--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if all(isnan(peak_firing_rate(:,unit))) || all(peak_firing_rate(:,unit) < 1)
        continue %unit doesn't fire enough go to next unit
    elseif all((spatial_info.shuffled_rate_prctile(:,unit) < 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) < 95))
        continue %unit likely not spatial
    elseif  ~any((spatial_info.shuffled_rate_prctile(:,unit) > 95) & (spatial_info.shuffled_spatialstability_prctile(:,unit) > 95))
        continue %only process cells that pass both criterion under at least 1 condition
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate View Field Location---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %draw spikes in locations above threshold (currently > 20% of max)
    which_condition_skaggs = spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95;
    which_condition_stable = spatial_info.shuffled_spatialstability_prctile(:,unit) >= 95;
    which_condition_both = which_condition_skaggs & which_condition_stable;
    
    if which_condition_both(3) ~= 1
        continue
    end
    
    filtered_time = filter_time(eyepos{unit},imageX,imageY,Fs,binsize,H);
    filtered_time(filtered_time < min_bin_dur) = NaN; %can cause aribitrarily high firing rates 25+ ms or more
    filtered_space = filter_space(eyepos{unit},spike_times{unit},imageX,imageY,binsize,H);
    
    firing_rate = filtered_space./filtered_time;
    fr = sort(firing_rate(1:end));
    fr(isnan(fr)) = [];
    maxfr = fr(round(0.95*length(fr)));% the ~95%-tile
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Define Place Field---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    sil = zeros(1,4); %determines the number of clusters by comparing the ratio
    %of intercluster and intracluster distances, faster mod of silhouette
    for numclusts = 2:4
        T = kmeans(fr',numclusts,'replicate',5);
        [silh] = InterVSIntraDist(fr',T);
        sil(numclusts) = mean(silh);
    end
    numclusters = find(sil == max(sil));
    T = kmeans(fr',numclusters(end),'replicate',25);
    
    region_mean_fr = zeros(1,numclusters);
    for t = 1:numclusters
        region_mean_fr(t) = mean(fr(T == t));
    end
    [~,place_field_cluster] = max(region_mean_fr);
    min_firing_rate = min(fr(T == place_field_cluster));
    
    [r,c] = find(firing_rate > min_firing_rate);
    firing_ind = sub2ind(size(firing_rate),r,c);
    threshold_matrix = zeros(size(firing_rate));
    threshold_matrix(firing_ind) = 1;
    threshold_matrix(isnan(firing_rate)) = NaN; %remove locations with too little data from further analysis
    if sum(sum(~isnan(threshold_matrix))) < 20
        continue %way too little area to process
    end
    threshold_matrix = imresize(threshold_matrix,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates
    area(unit) = 100*nansum(nansum(threshold_matrix))/sum(sum(~isnan(threshold_matrix)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Segment Fixations In/Out Field & Calculate Fixation Times---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %will need to calculate time/index in recording session for field trip analysis
    fixation_time = NaN(1,1000); %fixation start time within cfg structure starting from begining of session
    
    fixation_time_in_trial = NaN(1,1000); %fixation start time within trial
    trial_number = NaN(1,1000); %trial number for that fixation
    fix_in_out = NaN(1,1000); %firing rate for fixations inside (1) or outside (0) of view file
    
    fix_ind = 1; %fixation # so can track in variables above
    fixationstats = absolute_fixationstats; %reload because written over below
    cfg = absolute_cfg;
    num_trials = length(cfg.trl);
    %easier to start over rather than figure out format from imported data
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start+twin; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                imgoff = imgoff-twin; %to account for the twin when calculating firing rates and stuff
                
                %---Get eye data and remove fixations before/after image period---%
                fixationtimes = fixationstats{t}.fixationtimes;
                fixations = fixationstats{t}.fixations;
                %fixation started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(fixationtimes(1,:) > imgoff);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %remove fixations that are too short to really look at for analysis
                fixdur = fixationtimes(2,:)-fixationtimes(1,:)+1;
                fixationtimes(:,fixdur < minimum_fix_duration) = [];
                fixations(:,fixdur < minimum_fix_duration) = [];
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                if length(img_index) > 1
                    imgnum = which_img(img_index(1));
                    img_index = find(which_img == imgnum);
                    which_img(img_index) = NaN;
                    img_cnd(img_index) = NaN;
                    novel_vs_repeat(img_index) = NaN;
                    img_index = find(img_cnd == cfg.trl(t).cnd);
                end
                if isempty(img_index)
                    continue
                end
                
                for f = 1:size(fixations,2)
                    if f > 1
                        %determine if prior fixation was in our out of view field
                        fixx = round(fixations(1,f-1));
                        fixx(fixx < 1) = 1;
                        fixx(fixx > imageX) = imageX;
                        fixy = imageY-round(fixations(2,f-1));
                        fixy(fixy < 1) = 1;
                        fixy(fixy > imageY) = imageY;
                        
                        if threshold_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                            fixx = round(fixations(1,f));
                            fixx(fixx < 1) = 1;
                            fixx(fixx > imageX) = imageX;
                            fixy = imageY-round(fixations(2,f));
                            fixy(fixy < 1) = 1;
                            fixy(fixy > imageY) = imageY;
                            if threshold_matrix(fixy,fixx) == 1 %then inside
                                fix_in_out(fix_ind) = 2;%want to keep for spike-field coherence analysis
                                fixt = fixationtimes(1,f);
                                fixation_time(fix_ind) = fixt+trial_start; %start of fixation within trial
                                fixation_time_in_trial(fix_ind) = fixt;% start of fixation in file
                                trial_number(fix_ind) = t; %trial number for that fixation
                                fix_ind = fix_ind+1;
                                continue %move to the next one
                            end
                        end
                    end
                    
                    %determine if fixation was in our out of view field
                    fixx = round(fixations(1,f));
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-round(fixations(2,f));
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if threshold_matrix(fixy,fixx) == 1 %then inside
                        fix_in_out(fix_ind) = 1;
                    elseif threshold_matrix(fixy,fixx) == 0 %then inside, NaNs are for locations not analyzed
                        fix_in_out(fix_ind) = 0;
                    else %not a valid fixation location too sparse of coverage to analyze
                        continue
                    end
                    
                    fixt = fixationtimes(1,f);
                    fixation_time(fix_ind) = fixt+trial_start; %start of fixation within trial
                    fixation_time_in_trial(fix_ind) = fixt;% start of fixation in file
                    trial_number(fix_ind) = t; %trial number for that fixation
                    
                    fix_ind = fix_ind+1;
                end
            end
        end
    end
    
    %remove excess NaNs associated with error trials
    fix_in_out = laundry(fix_in_out);
    fixation_time = laundry(fixation_time);
    trial_number = laundry(trial_number);
    fixation_time_in_trial = laundry(fixation_time_in_trial);
    
    
    %---Calculate Fixation Locked Firing and STA LFP---%
    fix_locked_firing =  NaN(length(fix_in_out),2*twin); %firing rate locked to fixations
    STA_LFP = cell(1,4); %need a cell for each channel
    STA_LFP(:) = {zeros(length(fix_in_out),2*twin)};
    
    %code to find spike/LFP channels
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
    goodchan = find(~isnan(LFPchannels));
    goodchan = goodchan(1);
    
    for f = 1:length(fix_in_out);
        %for firing rate
        spikes = find(data(unit).values{trial_number(f)});
        fixt = fixation_time_in_trial(f);
        fix_spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin)-fixt+twin;
        temp = zeros(1,twin*2);
        temp(fix_spikes) = 1;
        fix_locked_firing(f,:) = temp;
        
        %for STA
        if isempty(fix_spikes)
            continue
        end
        
        spikes = spikes(spikes > fixt-twin & spikes <= fixt+twin);
        %remove spikes that wont have LFP data also image probably won't be
        %on by then anyway
        spikes(spikes+twin > length(data(LFPchannels(goodchan)).values{trial_number(f)})) = [];
        if isempty(spikes)
            for l = 1:length(LFPchannels)
                for s = 1:length(spikes);
                    STA_LFP{l}(f,:) = NaN; %don't count if didn't fire for that fixation
                end
            end
        else
            total_spikes = length(spikes);
            for l = 1:length(LFPchannels)
                if ~isnan(LFPchannels(l))
                    LFP = data(LFPchannels(l)).values{trial_number(f)};
                    for s = 1:length(spikes);
                        STA_LFP{l}(f,:) = STA_LFP{l}(f,:)+LFP(spikes(s)-twin+1:spikes(s)+twin)/total_spikes;
                    end
                end
            end
        end
    end
    
    LFPchannel_nums = 1:4;
    LFPchannel_nums(isnan(LFPchannels)) = [];
    LFPchannels(isnan(LFPchannels)) = []; %remove now for fieldtrip was conserving for STA format
    
    %store variables across units for later access
    list_fixation_locked_firing{unit} = fix_locked_firing;
    in_out{unit} = fix_in_out;
    list_STA_LFP{unit} = STA_LFP;
   

    %%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot---%%%
    %%%%%%%%%%%%%%%%%%%%%
    figure
    
    %draw raw spikes on all images all locations
    make_spike_jittered_plot(eyepos{unit},spike_times{unit},[4 4],1)
    title('Raw: All images')
    set(gca,'Xcolor','w')
    set(gca,'Ycolor','w')
    
    %draw firing rate map
    subplot(4,4,2)
    h = imagesc(firing_rate);
    set(h,'alphadata',~isnan(filtered_time));
    colormap('jet');
    title('Rate Map: All images')
    axis off
    axis equal
    colorbar
    
    %---draw view field---%
    subplot(4,4,3)
    h = imagesc(threshold_matrix);
    axis off
    axis equal
    title(sprintf(['View Field Location \n Area = ' num2str(area(unit),2)]));
    
    %---Plot firing rate for in vs out of field---%
    t = -twin:twin-1;
    subplot(4,4,4)
    hold on
    dofill(t,fix_locked_firing(fix_in_out == 1,:),'blue',1,smval(1));%inside
    dofill(t,fix_locked_firing(fix_in_out == 0,:),'red',1,smval(1));%outside
    xlabel('Time from Fixation Onset (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Fix Inside','Fix Outside','Location','NorthWest')
    title(sprintf(['n_{in} = ' num2str(sum(fix_in_out == 1)) ', n_{out} = ' num2str(sum(fix_in_out == 0))]))
    
    %---Plot STA by channel and for in vs out of field---%
    for l = 1:length(STA_LFP)
        subplot(4,4,4+l)
        hold on
        plot(nanmean(STA_LFP{l}(fix_in_out >= 1,:)),'b') %inside
        plot(nanmean(STA_LFP{l}(fix_in_out == 0,:)),'r') %outside
        set(gca,'Xtick',0:twin/2:2*twin)
        set(gca,'XtickLabel',num2cell((0:twin/2:2*twin)-twin));
        ylabel('STA (a.u.)')
        xlabel('Time from Fixation Onset (ms)')
        
        if str2double(unit_names{unit}(7)) == l
            title(['Unit Channel: LFP chan ' num2str(l)])
        else
            title(['LFP chan ' num2str(l)])
        end
    end
    
    %---Plot Fourier Spectrum for Inside Field STA LFP---%
    Fs = 1000;            % Sampling frequency
    L = 2*twin;             % Length of signal
    n = 2^nextpow2(L);
    f = Fs*(0:(n/2))/n;
    
    
    for l = 1:length(STA_LFP)
        subplot(4,4,8+l)
        
        Y = fft(nanmean(STA_LFP{l}(fix_in_out >= 1,:)),n);%inside only
        P = abs(Y/n);
        plot(f,P(1:n/2+1))
        xlabel('f (Hz)')
        xlim([1 31])
        ylabel('Power')
    end
    
    %---Plot AutoCorrelations---%
    
    %for inside field
    fix_in = fix_locked_firing(fix_in_out >= 1,:);
    fix_in(sum(fix_in,2) == 0,:) = []; %remove trials without spikes
    in_xc = NaN(size(fix_in,1),2*size(fix_in,2)-1);
    for f = 1:size(fix_in,1)
        xc = xcorr(fix_in(f,:),fix_in(f,:));
        in_xc(f,:) = xc;
    end
    in_xc(:,twin*2) = 0; %remove 0 lag it will always be 1
    in_xc = in_xc/sum(sum(fix_in));
    
    %for outside field
    fix_out = fix_locked_firing(fix_in_out == 0,:);
    fix_out(sum(fix_out,2) == 0,:) = []; %remove trials without spikes
    out_xc = NaN(size(fix_out,1),2*size(fix_out,2)-1);
    for f = 1:size(fix_out,1)
        xc = xcorr(fix_out(f,:),fix_out(f,:));
        out_xc(f,:) = xc;
    end
    out_xc(:,twin*2) = 0; %remove 0 lag it will always be 1
    out_xc = out_xc/sum(sum(fix_out));
    
    t = (1:size(in_xc,2))-2*twin;
    subplot(4,4,13)
    hold on
    [~,~,~,y1] = dofill(t,out_xc,'red',1,smval(2));
    [~,~,~,y2] = dofill(t,in_xc,'blue',1,smval(2));
    xlim([-1000 1000])
    ylim([0 1.05*max([max(y1),max(y2)])])
    ylabel('AutoCorrelation (a.u.)')
    xlabel('lag')
    peak_lag =  abs(find(y1 == max(y1))-2*twin);
    peak_lag = peak_lag(1);
    title(['Peak: ' num2str(peak_lag) ' ms,' num2str(round(1000/peak_lag)) ' Hz'])
    
    subplot(4,4,14)
    hold on
    [~,~,~,y1] = dofill(t,out_xc,'red',1,smval(2));
    [~,~,~,y2] = dofill(t,in_xc,'blue',1,smval(2));
    xlim([-100 100])
    ylim([0 1.05*max([max(y1),max(y2)])])
    ylabel('AutoCorrelation (a.u.)')
    xlabel('lag')
    
    
    spike_autocorr{1,unit} = in_xc;
    spike_autocorr{2,unit} = out_xc;
    
    %---Plot Fourier Spectrum for Inside Field AutoCorr---%
    Fs = 1000;            % Sampling frequency
    L = size(in_xc,2);             % Length of signal
    n = 2^nextpow2(L);
    f = Fs*(0:(n/2))/n;
    
    subplot(4,4,15)
    
    Y = fft(nanmean(in_xc),n);%inside only
    P = abs(Y/n);
    plot(f,P(1:n/2+1))
    xlabel('f (Hz)')
    xlim([1 21])
    ylabel('Power')
    
    subplot(4,4,16)
    P = abs(Y/n);
    plot(f,P(1:n/2+1))
    xlabel('f (Hz)')
    xlim([21 101])
    ylabel('Power')
    
    
    subtitle(['Spatial Plots' unit_names{unit}]);
    
    save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_Place_Cell_Freq_analysis']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Calculate Spike-Field Coherence---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    excludeChan       = str2num(unit_names{unit}(6)); % exclude the same channel
    chan              = true(1,length(LFPchannel_nums));
    chan(excludeChan==LFPchannel_nums) = false;
    param = 'ppc2'; % set the desired parameter, better for bursty low spike counts?
    
    
    % convert to field trip format [event start index-buffer1, event end index+buffer2, -buffer1]
    fix_in = fixation_time(fix_in_out >= 1); %takes any fixation in field even if not first
    fixation_aligned = [fix_in'-buffer fix_in'+buffer -buffer*ones(length(fix_in),1)];
    
    %remove trials if data extends past end of file
    fixation_aligned(fixation_aligned(:,2) > hdr.nSamples,:) = [];
    
    % read the data from file and preprocess them
    cfg.channel       = cfg.channel(LFPchannels)';
    cfg.dftfilter     = 'yes';
    cfg.dftfreq       = Fline;
    cfg.padding       = 1;
    cfg.continuous    = 'yes';
    
    %LFP data aligned saccades for list
    cfg.trl = fixation_aligned;
    fixation_aligneddata = ft_preprocessing(cfg);
    
    %read spike data into fieldtrip format
    spike = ft_read_spike(cfg.datafile);
    
    if all(strcmpi(spike.label,unit_names{unit}) == 0) %no spike found
       continue 
    end
    
    cfgs              = [];
    cfgs.spikechannel = {[unit_names{unit} '_wf']};
    spike            = ft_spike_select(cfgs, spike);
    
    %define trials based on cfg.trl above i.e. locked to fixations
    cfgs           = [];
    cfgs.hdr       = fixation_aligneddata.hdr; % contains information for conversion of samples to timestamps
    cfgs.trlunit   = 'samples';
    cfgs.trl       = cfg.trl; % now in samples
    spikeTrials   = ft_spike_maketrials(cfgs,spike);
    
    if size(spikeTrials.timestamp{1},2) <= 1%something odd with spike read
       continue 
    end
   

    %---Calculate for Low Frequencies---%
    cfg           = [];
    cfg.method    = 'mtmconvol';
    cfg.foi         = 2:1:30;
    cfg.t_ftimwin   = 5./cfg.foi;  % 7 cycles per time window
    cfg.taper     = 'hanning';
    stsConvol    = ft_spiketriggeredspectrum(cfg, fixation_aligneddata, spikeTrials);

    k = 1;
    % compute the statistics on the phases
    cfg               = [];
    cfg.method        = param; % compute the Pairwise Phase Consistency
    cfg.spikechannel  = stsConvol.label{k};
    cfg.channel       = stsConvol.lfplabel(chan); % selected LFP channels
    cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    cfg.timwin        = 'all'; % compute over all available spikes in the window
    cfg.latency       = [-0.5 0.5]; % sustained visual stimulation period
    low_statSts       = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
    all_low_statSts{unit} = low_statSts; 
    
    % plot the results
    figure
    subplot(2,2,1)
    plot(low_statSts.freq,low_statSts.ppc2')
    xlim([2 30])
    xlabel('frequency')
    ylabel('PPC')
    title('Low Frequency')
    
%     cfg                = [];
%     cfg.method         = param;
%     cfg.spikechannel   = stsConvol.label{k};
%     cfg.channel        = stsConvol.lfplabel(chan);
%     cfg.avgoverchan    = 'unweighted';
%     cfg.winstepsize    = 0.01; % step size of the window that we slide over time
%     cfg.timwin         = 0.5; % duration of sliding window
%     low_statSts            = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
%     
%     low_statSts.(param) = permute(conv2(squeeze(low_statSts.(param)), ones(1,20)./20, 'same'),[3 1 2]); % apply some smoothing over 0.2 sec
%     
%     cfg            = [];
%     cfg.parameter  = param;
%     cfg.refchannel = low_statSts.labelcmb{1,1};
%     cfg.channel    = low_statSts.labelcmb{1,2};
%     cfg.xlim       = [-1 2];
%     
%     subplot(2,2,3)
%     imagesc(low_statSts.time,low_statSts.freq,squeeze(low_statSts.ppc2))
%     hold on
%     plot([0 0],[low_statSts.freq(1) low_statSts.freq(end)],'w--')
%     hold off
%     axis xy
%     ylabel('Frequency (Hz)')
%     xlim([-0.5 0.5])
%     xlabel('Time from Fixation (sec)')

    
    
    %---Calculate for high Frequencies---%
    cfg           = [];
    cfg.method    = 'mtmconvol';
    cfg.foi         = 30:3:90;
    cfg.t_ftimwin   = 5./cfg.foi;  % 7 cycles per time window
    cfg.taper     = 'hanning';
    stsConvol    = ft_spiketriggeredspectrum(cfg, fixation_aligneddata, spikeTrials);
    
    k = 1;
    % compute the statistics on the phases
    cfg               = [];
    cfg.method        = param; % compute the Pairwise Phase Consistency
    cfg.spikechannel  = stsConvol.label{k};
    cfg.channel       = stsConvol.lfplabel(chan); % selected LFP channels
    cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    cfg.timwin        = 'all'; % compute over all available spikes in the window
    cfg.latency       = [-0.5 0.5]; % sustained visual stimulation period
    high_statSts      = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
    all_high_statSts = high_statSts;
    
    % plot the results
    subplot(2,2,2)
    plot(high_statSts.freq,high_statSts.ppc2)
    xlabel('frequency')
    xlim([30 90])
    ylabel('PPC')
    title('High Frequency')
    
%     param = 'ppc2'; % set the desired parameter
%     cfg                = [];
%     cfg.method         = param;
%     cfg.spikechannel   = stsConvol.label{k};
%     cfg.channel        = stsConvol.lfplabel(chan);
%     cfg.avgoverchan    = 'unweighted';
%     cfg.winstepsize    = 0.01; % step size of the window that we slide over time
%     cfg.timwin         = 0.5; % duration of sliding window
%     high_statSts            = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
%     
%     high_statSts.(param) = permute(conv2(squeeze(high_statSts.(param)), ones(1,20)./20, 'same'),[3 1 2]); % apply some smoothing over 0.2 sec
%     
%     cfg            = [];
%     cfg.parameter  = param;
%     cfg.refchannel = high_statSts.labelcmb{1,1};
%     cfg.channel    = high_statSts.labelcmb{1,2};
%     cfg.xlim       = [-1 2];
%     
%     subplot(2,2,4)
%     imagesc(high_statSts.time,high_statSts.freq,squeeze(high_statSts.ppc2))
%     hold on
%     plot([0 0],[high_statSts.freq(1) high_statSts.freq(end)],'w--')
%     hold off
%     axis xy
%     ylabel('Frequency (Hz)')
%     xlim([-0.5 0.5])
%     xlabel('Time from Fixation (sec)')
    
    save_and_close_fig(view_cell_dir,[task_file(1:end-11) '_' unit_names{unit} '_Place_Cell_SpikeField_analysis']);
    
end
save([data_dir task_file(1:8) '-Place_Cell_Freq.mat'],...
    'twin','smval','list_fixation_locked_firing','in_out',...
    'task_file','area','list_STA_LFP','spike_autocorr',...
    'all_low_statSts','all_high_statSts')
end