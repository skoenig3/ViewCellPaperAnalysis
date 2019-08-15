function ListSQ_Saccade_Direction_and_Amplitude_Analysis_Future_and_Past(data_dir,figure_dir,session_data)
%written by Seth Konig 4/6/17
%written as a combination of code from ListSQ_Saccade_Amplitude_Analysis,
%ListSQ_Saccade_Direction_Analysis, and ListSQ_Saccade_UpCommingDirection_Analysis.
%
%It became aparent that saccade alignment was a better predictor of these
%than fixation alignment. Further, it became clear that some neurons encode
%the upcomming saccade direction, some encode peri-saccadic direction
%infomration (including during), and some encode post-saccadic direction
%information. Since redoing anlaysis I've also added amplitude to this
%analysis as well since it could be equally effected.

%1) code determines whether units are modulation by saccade direction. Currently
%For units that are spatially modulated (95% corr 1/2 &  95% skaggs)
% the data is split into in field and out of field regions. Code calculates
% significance with circ_rtest and by bootstrapping (permutation) of mrl. mrl
% seems more accurate in the sense that it taskes into account sample size
% better and should be more sensitive to low firing rate neurons.
%
%2)%Code determines whether units are modulation by saccade ampltiude. Same
%as above!

% code rechecked for bugs SDK 5/2/17

figure_dir = [figure_dir 'Saccade Direction and Amplitude\'];
task = 'ListSQ';

%---generic parameters---%
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
image_on_twin = 500;%don't use data within 500 ms of fixation onset due to strong visual response
min_window_width = 50; %time around peak to take for firing rate, based on FWHM of place cell fixation analysis
numshuffs = 10000; %number of times to shuffled/permutate saccade direction analysis
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
ITIstart_code = 15; %start of ITI/trial
img_on_code = 23; %image turned on
img_off_code = 24; %image turned off
smval = 30 ;%15 ms std, want to smooth at high frequency modulations
min_fix_dur = 100; %100 ms, don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva, don't want mini/micro saccades too small and hard to detect
min_saccades_with_spikes = 0.05;% 5% to remove neurons that have basically no activity since they pass too :(...
%essentially 0.4 Hz threshold if window is 50 ms in wide

%---saccade direction and amplitude parameters---%
bin_deg = 3;%1 %number of degrees per bin
bin_deg2 = 45; %initial degree bins to determine time of peak direction modualtion
smval_deg = 6;%18 %9 degrees std
bin_amplitude = 2; %dva for calculating saccade amplitude tuning
bin_amplitude2 = 4;%dva for estimating window of interest
max_amplitude = 16;%dva, approximately 95th percentile, don't have many large ones so remove

twinad1 = 200;%presaccade window
twinad2 = 400;%pos saccade window
median_sac_dur = 44;
tm = -twinad1+1:twinad2;
start_window = twinad1-min_fix_dur;%how early before saccade can you look for direction tuning
end_window = twinad1+min_fix_dur+median_sac_dur;%how late after saccade can you look for direction tuning,
%some neurons have large latencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Load important Session Data and Information---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = 'ListSQ';
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    disp('No ListSQ file could be found. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr','fixationstats');

%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return %since no units skip analysis
end

%---determine number of blocks unit was stable for---%
%remove units with too few trials
%these are the absolute minimum data required to do data analysis
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end
disp([task_file ': running Saccade modulation anlaysis...'])

load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
    'spatial_info','eyepos','spike_times','binsize','Fs','filter_width')
H = define_spatial_filter(filter_width);
load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],'all_place_field_matrix')


%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
Fs = data(1).fsample; %should be 1000

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,~,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
num_trials = length(cfg.trl);%number of trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Get Saccade Eye Movement Information---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saccade_aligned_firing = cell(1,num_units);
saccade_directions = cell(1,num_units);
saccade_amplitudes = cell(1,num_units);
fixation_starts = cell(1,num_units);
fixation_ends = cell(1,num_units);
sac_in_out = cell(1,num_units);
%1) first fixation in: out-> in
%2) fixation in field but not first: in -> in
%3) first fixation out of field: in -> out
%4) fixation out of field but not first: out-> out


for unit = 1:num_units
    if isempty(eyepos{unit})
        continue %no data for this neuron
    end
    
    %---Calculate Firing Rate Locked to Saccades---%
    sac_in_out{unit} = NaN(1,3000); %firing rate for fixations inside or outside  of place field
    saccade_aligned_firing{unit} = NaN(3000,(twinad1+twinad2)); %spike trains locked to saccade start
    saccade_directions{unit} = NaN(1,3000);%saccade directions
    saccade_amplitudes{unit} = NaN(1,3000);%saccade amplitudes
    fixation_starts{unit} = NaN(1,3000);%start of prior fixation relative to saccade start
    fixation_ends{unit} = NaN(1,3000); %end of following fixation relative to saccade start
    
    sac_ind = 1; %saccade # so can track in variables above
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                
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
                invalid= find(fixationtimes(1,:) > imgoff-twinad2);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %saccade started before image turned on
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(saccadetimes(1,:) > imgoff-twinad2);
                saccadetimes(:,invalid) = [];
                
                %remove fixations that are too short to really look at for analysis
                fixdur = fixationtimes(2,:)-fixationtimes(1,:)+1;
                
                img_index = find(img_cnd == cfg.trl(t).cnd);
                %which image monkey viewed if nan or empty then skip cuz
                %bad trial
                if isempty(img_index) || any(isnan(which_img(img_index)))
                    continue
                end
                
                spikes = find(data(unit).values{t}); %spike trains for this trial
                for s = 1:size(saccadetimes,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                    prior_fix = find(saccadetimes(1,s) == fixationtimes(2,:)+1);%find preceding fixation
                    post_fix =  find(saccadetimes(2,s) == fixationtimes(1,:)-1);%find following fixation
                    
                    if isempty(prior_fix) || isempty(post_fix) %no prior or following fixation so was proabbly looking off screen
                        continue;
                    end
                    
                    pre_dur = fixdur(prior_fix);
                    post_dur = fixdur(post_fix);
                    
                    if pre_dur < min_fix_dur || post_dur < min_fix_dur
                        continue
                    end
                    
                    sacamp = sqrt(sum((xy(:,saccadetimes(2,s))-xy(:,saccadetimes(1,s))).^2)); %saccade amplitude
                    if sacamp < min_sac_amp %prior saccade is too small so ignore
                        continue
                    end
                    saccade_amplitudes{unit}(sac_ind) = sacamp;
                    
                    prior_fix_in_out = NaN;
                    %determine if prior fixation was in or out of place field
                    fixx = fixations(1,prior_fix);
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-fixations(2,prior_fix);
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if all_place_field_matrix{unit}(fixy,fixx) == 1 %then prior fixation was already inside
                        prior_fix_in_out = 1;
                    else
                        prior_fix_in_out = 0;
                    end
                    last_fixx = fixx;
                    last_fixy = fixy;
                    
                    %determine if this fixation was in our out of place field
                    fixx = fixations(1,post_fix);
                    fixx(fixx < 1) = 1;
                    fixx(fixx > imageX) = imageX;
                    fixy = imageY-fixations(2,post_fix);
                    fixy(fixy < 1) = 1;
                    fixy(fixy > imageY) = imageY;
                    if all_place_field_matrix{unit}(fixy,fixx) == 1 %then inside
                        if prior_fix_in_out == 1%prior fixation was inside so in->in
                            sac_in_out{unit}(sac_ind) = 2;
                        else %out->in
                            sac_in_out{unit}(sac_ind) = 1;
                        end
                    elseif all_place_field_matrix{unit}(fixy,fixx) == 0 %then inside
                        if prior_fix_in_out == 1%prior fixation was inside so in->out
                            sac_in_out{unit}(sac_ind) = 3;
                        else %out->out
                            sac_in_out{unit}(sac_ind) = 4;
                        end
                    else %not a valid fixation location too sparse of coverage to analyze,
                        %NaNs are for locations not analyzed
                        continue
                    end
                    
                    
                    %get firing rate locked to fixation
                    sact = saccadetimes(1,s);%start of fixation
                    sac_spikes = spikes(spikes > sact-twinad1 & spikes <= sact+twinad2)-sact+twinad1;
                    temp = zeros(1,twinad1+twinad2);
                    temp(sac_spikes) = 1;
                    saccade_aligned_firing{unit}(sac_ind,:) = temp;
                    
                    %reflip y otherwise direction gets flipped too
                    fixy = fixations(2,post_fix);
                    last_fixy = fixations(2,prior_fix);
                    saccade_directions{unit}(sac_ind) = atan2d(fixy-last_fixy,fixx-last_fixx);
                    
                    fixation_starts{unit} (sac_ind) =fixationtimes(1,prior_fix)-sact;
                    fixation_ends{unit}(sac_ind) = fixationtimes(2,post_fix)-sact;
                    
                    sac_ind = sac_ind+1;
                end
            end
        end
    end
end
%%
%---remove excessive NaNs---%
saccade_aligned_firing = laundry(saccade_aligned_firing);
saccade_directions = laundry(saccade_directions);
saccade_amplitudes = laundry(saccade_amplitudes);
sac_in_out = laundry(sac_in_out);
fixation_starts = laundry(fixation_starts);
fixation_ends = laundry(fixation_ends);


%%
%---Where to Store Direction Summary Data---%
mrls.all_saccades = NaN(1,num_units); %MRLs all fixations ignoring out2in and in2out
mrls.all_saccades_shuffled = cell(1,num_units); %shuffled MRLs all fixations ignoring out2in and in2out
mrls.all_saccades_shuffled_prctile = NaN(1,num_units); %observed MRL percentile
mrls.in2in = NaN(1,num_units);  %MRLs for in2in fixations ignoring all others
mrls.in2in_shuffled = cell(1,num_units);  %shuffled MRLs for in2in fixations ignoring all others
mrls.in2in_shuffled_prctile = NaN(1,num_units);  %observed MRLs percentile for in2in fixations ignoring all others
mrls.out2out = NaN(1,num_units);  %MRLs for out2out fixations ignoring all others
mrls.out2out_shuffled = cell(1,num_units); %shuffled MRLs for out2out fixations ignoring all others
mrls.out2out_shuffled_prctile = NaN(1,num_units);%MRLs percentile for out2out fixations ignoring all others
uniformity_pvalue = NaN(3,num_units); %circ_rtest...row 1: all fixations, row 2: in2in, row3: out2out
direction_binned_firing_rate_curves = cell(3,num_units); %row 1: all fixations, row 2: in2in, row3: out2out
all_direction_windows = cell(1,num_units);%window used to calculate MRLs

%---Where to Store Amplitude Summary Data---%
amplitude_correlations = NaN(1,num_units);%observed correlation between firing rate and saccade amplitude
shuffled_amplitude_correlations = cell(1,num_units); %shuffled correlations
amplitude_correlations_percentile = NaN(1,num_units); %observed correlation percentile
amplitude_binned_firing_rate_curves = cell(1,num_units); %"firing rate curves"
all_amplitude_windows = cell(1,num_units);%window used to calculate amplitude tuning


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Determine if Units are Modulation by Direction or Amplitude---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unit = 1:num_units
    if isempty(saccade_aligned_firing{unit}) %no valid trials so go to next unit
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Grab Saccade Aligned Activity--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sac_aligned = saccade_aligned_firing{unit}; %fixation aligned firing rate
    sac_dirs = saccade_directions{unit}; %saccade directions organized the same way as sac_algined
    sac_amps = saccade_amplitudes{unit}/24; %saccade amplitudes organized the same way as sac_algined
    fix_starts = fixation_starts{unit};
    fix_ends = fixation_ends{unit};
    
    %---Remove Counfounding Eye Movements for Spatially Modulated Neurons---%
    %remove saccades in2in or out2out since in2out or out2in
    %could be biased by field location creating artificial direction tuning
    %but only do this if spatial (both criterion)
    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        sac_aligned = sac_aligned(sac_in_out{unit} == 2 | sac_in_out{unit} == 4,:); %saccades in2in or out2out
        sac_dirs = sac_dirs(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for saccades in2in or out2out only
        sac_amps = sac_amps(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for saccades in2in or out2out only
        fix_starts = fix_starts((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
        fix_ends = fix_ends((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
    end
    
    %also remove saccades that are too big since won't have many anyway
    sac_aligned_amp = sac_aligned;
    too_large = find(sac_amps > max_amplitude);
    sac_aligned_amp(too_large,:) = [];
    sac_amps(too_large)=[];
    sac_fix_starts = fix_starts;
    sac_fix_starts(too_large) = [];
    sac_fix_ends = fix_ends;
    sac_fix_ends(too_large) = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Saccade Direction Analysis---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %---Crudely estimate direciton tuning---%
    fr = nandens(sac_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations
    degrees2 = [0:bin_deg2:360]-180;
    direction_curves = NaN(length(degrees2),twinad1+twinad2);
    for bin = 2:length(degrees2)
        these_directions = (sac_dirs < degrees2(bin) & sac_dirs >= degrees2(bin-1));
        direction_curves(bin,:) = nandens(sac_aligned(these_directions,:),smval,'gauss',Fs)-fr;
    end
    crude_direction_direction_curves = direction_curves; %save for later
    direction_curves = max(crude_direction_direction_curves)-min(crude_direction_direction_curves);%subtract lowest from highest firing rates
    direction_curves = nandens(direction_curves,round(smval/3),'gauss',1);%minor smoothing to remove any high frequency bumps
    
    %---get window above 1/2 max (Essentially FWHM)---%
    [max_tuning_ind,window] = get_half_max_window(direction_curves,start_window,end_window);
    if length(window) < min_window_width
        window = max_tuning_ind-min_window_width/2:max_tuning_ind+min_window_width/2;
    end
    
    %---Re-Estimate Direction Tuning More Accurately within Crude Window---%
    %will improve estimate as some neurons will have prefered tuning
    %directions in between the crude estimates above
    [mean_binned_firing_rate,degrees,~] = bin_directional_firing(bin_deg,sac_aligned,window,sac_dirs);
    [prefered_dirs,anti_prefered_dirs,~] = ...
        select_prefred_indeces(mean_binned_firing_rate,degrees,sac_dirs,smval_deg,bin_deg);%estimate prefered and anti-prefered directions
    prefered_curve = nandens(sac_aligned(prefered_dirs,:),smval,'gauss',Fs);
    anti_prefered_curve = nandens(sac_aligned(anti_prefered_dirs,:),smval,'gauss',Fs);
    pref_anti_curve = prefered_curve-anti_prefered_curve;
    [max_tuning_ind,window] = get_half_max_window(pref_anti_curve,start_window,end_window);
    if length(window) < min_window_width
        window = max_tuning_ind-min_window_width/2:max_tuning_ind+min_window_width/2;
    end
    all_direction_windows{unit} = window;
    
    %---Modify Window if Necessary based on Fixation Durations and Remove Fixation Durations that are too Short---%
    window_width = length(all_direction_windows{unit});
    window_end = all_direction_windows{unit}(end);
    window_start = all_direction_windows{unit}(1);
    if window_end > end_window
        window_end = window_end-twinad1;
        median_fixation_dur = median(fix_ends);
        %if more than half the fixations are shorter than window length
        %cut window down to median fixation duration othwerise we will just
        %cut fixations shorter than window end
        if median_fixation_dur < window_end
            window = all_direction_windows{unit};
            window(window > median_fixation_dur+twinad1) = [];
            if length(window) < min_window_width;
                all_direction_windows{unit} = [];
            else
                all_direction_windows{unit} = window;
            end
            window_end = all_direction_windows{unit}(end)-twinad1;
        end
        %remove fixations shorter than window as these could be
        %contaminated by the next saccade
        fixations_too_short = find(fix_ends < window_end);
        sac_dirs(fixations_too_short) = [];
        sac_aligned(fixations_too_short,:) = [];
        fix_ends(fixations_too_short) = [];
        fix_starts(fixations_too_short) = [];
    end
    if window_start < twinad1-min_fix_dur
        fixations_too_short = find((twinad1+fix_starts) > window_start);
        sac_dirs(fixations_too_short) = [];
        sac_aligned(fixations_too_short,:) = [];
        fix_ends(fixations_too_short) = [];
        fix_starts(fixations_too_short) = [];
    end
    
    if ~isempty(all_direction_windows{unit})
        
        %---Calculate Observed Saccade Direction Tuning for All Fixations---%
        window_width = length(all_direction_windows{unit});
        [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_aligned,all_direction_windows{unit},sac_dirs);
        mrls.all_saccades(unit) = mrl; %observed MRL
        direction_binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
        uniformity_pvalue(1,unit) = circ_rtest(sac_dirs*pi/180,sum(sac_aligned(:,all_direction_windows{unit}),2)*1000/window_width);
        %circ p_values may not be that useful since depends on firing rate and sample size
        
        spike_counts = sum(sac_aligned(:,all_direction_windows{unit}),2);
        if sum(spike_counts > 0) > min_saccades_with_spikes*size(sac_aligned,1)
            
            %---Calculate Bootstrapped MRLs for All Fixations---%
            shuffled_mrls = NaN(1,numshuffs); %shuffled mrls
            this_window = all_direction_windows{unit};
            spike_counts = sum(sac_aligned(:,this_window),2);
            parfor shuff = 1:numshuffs
                ind = randperm(length(sac_dirs)); %shuffle indeces
                dirs = sac_dirs(ind); %shuffled saccade directions
                mrl = circ_r(dirs*pi/180,spike_counts); %MRL for unbinned data
                shuffled_mrls(shuff) = mrl;
            end
            mrls.all_saccades_shuffled{unit} = shuffled_mrls;
            mrls.all_saccades_shuffled_prctile(unit) = 100*sum(mrls.all_saccades(unit) > shuffled_mrls)/numshuffs; %shuffled percentile
            
            %---Saccade Direction across Fixations In2In or Out2Out---%
            %only run if spatial (both criterion) or possibly spatial (either criterion)
            if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                
                %---In2In Fixations---%
                sac_in_in = saccade_aligned_firing{unit}(sac_in_out{unit} == 2,:); %fixation aligned firing rate
                dir_in_in = saccade_directions{unit}(sac_in_out{unit} == 2); %saccade directions
                fix_starts = fixation_starts{unit}(sac_in_out{unit} == 2);
                fix_ends = fixation_ends{unit}(sac_in_out{unit} == 2);
                
                %---Remove Fixations with Durations that are too short---%
                window_end = all_direction_windows{unit}(end);
                window_start = all_direction_windows{unit}(1);
                if window_end > end_window
                    fixations_too_short = find(fix_ends < window_end-twinad1);
                    dir_in_in(fixations_too_short) = [];
                    sac_in_in(fixations_too_short,:) = [];
                    fix_ends(fixations_too_short) = [];
                    fix_starts(fixations_too_short) = [];
                end
                if window_start < twinad1-min_fix_dur
                    fixations_too_short = find((twinad1+fix_starts) > window_start);
                    dir_in_in(fixations_too_short) = [];
                    sac_in_in(fixations_too_short,:) = [];
                    fix_ends(fixations_too_short) = [];
                    fix_starts(fixations_too_short) = [];
                end
                
                %---Calculate Observed Saccade Direction Tuning for In2In Fixations---%
                [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_in_in,all_direction_windows{unit},dir_in_in);
                mrls.in2in(unit) = mrl; %observed MRL
                direction_binned_firing_rate_curves{2,unit} = mean_binned_firing_rate; %binned firing rate
                uniformity_pvalue(2,unit) = circ_rtest(dir_in_in*pi/180,sum(sac_in_in(:,all_direction_windows{unit}),2)*1000/window_width);
                %circ p_values may not be that useful since depends on firing rate and
                %sample size
                
                %---Calculate Bootstrapped MRLs for In2In Fixations---%
                shuffled_mrls = NaN(1,numshuffs);
                spike_counts = sum(sac_in_in(:,this_window),2);
                parfor shuff = 1:numshuffs
                    ind = randperm(length(dir_in_in)); %shuffle indeces
                    dirs = dir_in_in(ind); %shuffled directions
                    mrl = circ_r(dirs*pi/180,spike_counts); %MRL for unbinned data
                    shuffled_mrls(shuff) = mrl;
                end
                mrls.in2in_shuffled{unit} = shuffled_mrls;
                mrls.in2in_shuffled_prctile(unit) = 100*sum(mrls.in2in(unit) > shuffled_mrls)/numshuffs; %shuffled percentile
                
                %---Out2Out Fixations---%
                sac_out_out = saccade_aligned_firing{unit}(sac_in_out{unit} == 4,:); %fixation aligned firing rate
                dir_out_out = saccade_directions{unit}(sac_in_out{unit} == 4); %saccade directions
                fix_starts = fixation_starts{unit}(sac_in_out{unit} == 4);
                fix_ends = fixation_ends{unit}(sac_in_out{unit} == 4);
                
                %---Remove Fixations with Durations that are too short---%
                window_end = all_direction_windows{unit}(end);
                window_start = all_direction_windows{unit}(1);
                if window_end > end_window
                    fixations_too_short = find(fix_ends < (window_end-twinad1));
                    dir_out_out(fixations_too_short) = [];
                    sac_out_out(fixations_too_short,:) = [];
                    fix_ends(fixations_too_short) = [];
                    fix_starts(fixations_too_short) = [];
                end
                if window_start < twinad1-min_fix_dur
                    fixations_too_short = find((twinad1+fix_starts) > window_start);
                    dir_out_out(fixations_too_short) = [];
                    sac_out_out(fixations_too_short,:) = [];
                    fix_ends(fixations_too_short) = [];
                    fix_starts(fixations_too_short) = [];
                end
                
                %---Calculate Observed Saccade Direction Tuning for Out2Out Fixations---%
                [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_out_out,all_direction_windows{unit},dir_out_out);
                mrls.out2out(unit) = mrl; %observed MRL
                direction_binned_firing_rate_curves{3,unit} = mean_binned_firing_rate; %binned firing rate
                uniformity_pvalue(3,unit) = circ_rtest(dir_out_out*pi/180,sum(sac_out_out(:,all_direction_windows{unit}),2)*1000/window_width);
                %circ p_values may not be that useful since depends on firing rate and
                %sample size
                
                %---Calculate Bootstrapped MRLs for Out2Out Fixations---%
                shuffled_mrls = NaN(1,numshuffs);
                spike_counts = sum(sac_out_out(:,this_window),2);
                parfor shuff = 1:numshuffs
                    ind = randperm(length(dir_out_out)); %shuffle indeces
                    dirs = dir_out_out(ind);%shuffled saccade directions
                    mrl = circ_r(dirs*pi/180,spike_counts); %MRL for unbinned data
                    shuffled_mrls(shuff) = mrl;
                end
                mrls.out2out_shuffled{unit} = shuffled_mrls;
                mrls.out2out_shuffled_prctile(unit) = 100*sum(mrls.out2out(unit) > shuffled_mrls)/numshuffs;%shuffled perecentile
            end
        else
            window_width = length(all_direction_windows{unit});
            [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_aligned,all_direction_windows{unit},sac_dirs);
            mrls.all_saccades(unit) = mrl; %observed MRL
            direction_binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
            uniformity_pvalue(1,unit) = circ_rtest(sac_dirs*pi/180,sum(sac_aligned(:,all_direction_windows{unit}),2)*1000/window_width);
            %circ p_values may not be that useful since depends on firing rate and sample size
            disp('Too sparse')
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Saccade Direction Plots---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    %---Fixation Aligned Activity out2out & in2in---%
    subplot(3,3,1)
    dofill(tm,sac_aligned,'black',1,smval); %smoothed fiirng rate curve
    hold on
    yl = ylim;
    if yl(1) < 0
        ylim([0 yl(2)]);
        yl(1) = 0;
    end
    if ~isempty(all_direction_windows{unit})
        h = fill([all_direction_windows{unit}(1) all_direction_windows{unit}(end) ...
            all_direction_windows{unit}(end) all_direction_windows{unit}(1) all_direction_windows{unit}(1)]-twinad1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
        set(h,'facealpha',.25,'EdgeColor','None')
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-twinad1 twinad2]);
    xlabel('Time from Saccade Start (ms)')
    ylabel('Firing Rate (Hz)')
    title('All Saccade Aligned Activity')
    
    %---Plot Firing Rate Curves for each of the cardinal and intermediate directions---%
    subplot(3,3,4)
    plot(tm,crude_direction_direction_curves')
    hold on
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-twinad1 twinad2]);
    xlabel('Time from Saccade Start (ms)')
    ylabel('Firing Rate (Hz)')
    title(['Firing Rate Curves by Crude 45' char(176) ' Direction']);
    box off
    
    %---Fixation Aligned Raster Sorted by Saccade Direction for out2out & in2in fixations---%
    subplot(3,3,5)
    [~,si] = sort(sac_dirs);
    sac_aligned_sorted = sac_aligned(si,:);
    [trial,time] = find(sac_aligned_sorted == 1);
    plot(time-twinad1,trial,'.k')
    xlim([-twinad1 twinad2])
    if ~isempty(trial)
        yl = [0 max(trial)+1];
        ylim([0 max(trial)+1]);
        hold on
        if ~isempty(all_direction_windows{unit})
            h = fill([all_direction_windows{unit}(1) all_direction_windows{unit}(end) ...
                all_direction_windows{unit}(end) all_direction_windows{unit}(1) all_direction_windows{unit}(1)]-twinad1,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
            set(h,'facealpha',.25,'EdgeColor','None')
        end
        hold off
    end
    xlim([-twinad1 twinad2]);
    ylabel('Ranked Sac. Dir.')
    xlabel('Time from Saccade Start (ms)')
    title('Fixation Aligned-Saccade Direction')
    box off
    
    %---Firing Rate Map---%
    subplot(3,3,6)
    firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
    maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
    h = imagesc(firing_rate_map);
    set(h,'alphadata',~isnan(firing_rate_map));
    axis off
    axis equal
    colormap('jet')
    colorbar
    clim = caxis;
    caxis([clim(1) maxfr])
    title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(unit) > 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
    end
    if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
    end
    if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
        title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
    if ~isempty(all_direction_windows{unit})
        %---Fixation Aligned Acitivity out2out & in2in Plotted by Prefered and AntiPrefered Direction-All Fixations---%
        if ~isempty(all_direction_windows{unit})
            [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
                select_prefred_indeces(direction_binned_firing_rate_curves{1,unit},degrees,sac_dirs,smval_deg,bin_deg);
            subplot(3,3,2)
            hold on
            dofill(tm,sac_aligned(prefered_dirs,:),'green',1,smval); %plot fixation aligned activity for prefered direction
            dofill(tm,sac_aligned(anti_prefered_dirs,:),'black',1,smval); %plot fixation aligned activity for anti-prefered direction
            yl = ylim;
            if yl(1) < 0;
                yl(1) = 0;
                ylim([0 yl(2)]);
            end
            plot([0 0],[yl(1) yl(2)],'k--')
            hold off
            xlim([-twinad1 twinad2]);
            xlabel('Time from Saccade Start (ms)')
            ylabel('Firing Rate (Hz)')
            title('All Saccade Aligned Activity')
            legend('Prefered','Anti-Preferred')
            
            %---Smoothed Polar Plot of Prefered Firing Direction-Fixations out2out & in2in---%
            subplot(3,3,3)
            polar(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b')
            if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                %only ran if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
                title(sprintf(['All mrl: ' num2str(mrls.all_saccades(unit),2) ' (' num2str(mrls.all_saccades_shuffled_prctile(unit),3) '%%)'...
                    ' \n In2In mrl: ' num2str(mrls.in2in(unit),2) ' (' num2str(mrls.in2in_shuffled_prctile(unit),3) '%%)']))
            else
                title(sprintf(['All mrl: ' num2str(mrls.all_saccades(unit),2) ' (' num2str(mrls.all_saccades_shuffled_prctile(unit),3) '%%)']))
            end
            
            %---Plots for Fixations Out2Out Only---%
            %not plotting In2In since this is less interesting and coverage isn't
            %great for all neurons anyway. I'm more intesteted in showing direction
            %tuning out of the field
            if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95) && ~isnan(mrls.out2out(unit))
                %only run if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
                
                %---Smoothed Plot of Firing Rate by Saccade Direction Out2Out ixations---%
                [prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
                    select_prefred_indeces(direction_binned_firing_rate_curves{3,unit},degrees,dir_out_out,smval_deg,bin_deg);
                subplot(3,3,7)
                polar(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b')
                title(['Out2Out mrl: ' num2str(mrls.out2out(unit),2) ' (' num2str(mrls.out2out_shuffled_prctile(unit),3) '%)'])
                
                %---Out2Out Fixation Aligned Raster Sorted by Saccade Direction---%
                subplot(3,3,8)
                [~,si] = sort(dir_out_out);
                sac_aligned_sorted = sac_out_out(si,:);
                [trial,time] = find(sac_aligned_sorted == 1);
                plot(time-twinad1,trial,'.k')
                xlim([-twinad1 twinad2])
                if ~isempty(trial)
                    yl = [0 max(trial)+1];
                    ylim([0 max(trial)+1]);
                    hold on
                    if ~isempty(all_direction_windows{unit})
                        h = fill([all_direction_windows{unit}(1) all_direction_windows{unit}(end) ...
                            all_direction_windows{unit}(end) all_direction_windows{unit}(1) all_direction_windows{unit}(1)]-twinad1,...
                            [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
                        set(h,'facealpha',.25,'EdgeColor','None')
                    end
                    hold off
                end
                xlim([-twinad1 twinad2]);
                ylabel('Ranked Sac. Dir.')
                xlabel('Time from Saccade Start (ms)')
                title('Out2Out Fixation Aligned-Saccade Direction')
                box off
                
                %---Fixation Aligned Acitivity Plotted by Prefered and AntiPrefered Direction-All Fixations---%
                subplot(3,3,9)
                hold on
                dofill(tm,sac_out_out(prefered_dirs,:),'green',1,smval);
                dofill(tm,sac_out_out(anti_prefered_dirs,:),'black',1,smval);
                yl = ylim;
                if yl(1) < 0;
                    yl(1) = 0;
                    ylim([0 yl(2)]);
                end
                plot([0 0],[yl(1) yl(2)],'k--')
                hold off
                xlim([-twinad1 twinad2]);
                xlabel('Time from Saccade Start (ms)')
                ylabel('Firing Rate (Hz)')
                title('Out2Out Fixation Aligned Activity')
                legend('Prefered','Anti-Preferred')%,'Neutral1','Neutral2')
            end
        end
        
        if multiunit(unit)
            subtitle(['Multiunit ' task_file(1:8) ' ' unit_stats{1,unit}]);
        else
            subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        end
        
    end
    
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Saccade_Direction_Analysis']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--Saccade Amplitude Anlaysis---%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amps = [0:bin_amplitude2:18]+2;
    amplitude_curves = NaN(length(amps),twinad1+twinad2);
    for bin = 1:length(amps)-1
        these_amplitudes = (sac_amps < amps(bin+1) & sac_amps >= amps(bin));
        amplitude_curves(bin,:) = nandens(sac_aligned_amp(these_amplitudes,:),smval,'gauss',Fs)-fr;
    end
    crude_amplitude_amplitude_curves = amplitude_curves; %save for later
    small = prctile(sac_amps,20);
    small = nandens(sac_aligned_amp(sac_amps <= small,:),smval,'gauss',Fs);
    large = prctile(sac_amps,80);
    large = nandens(sac_aligned_amp(sac_amps >= large,:),smval,'gauss',Fs);
    amp_diff = large-small;
    if abs(min(amp_diff(start_window:end_window))) > max(amp_diff(start_window:end_window))
        [min_tuning,max_tuning_ind] = min(amp_diff(start_window:end_window));
        max_tuning_ind = max_tuning_ind+start_window-1;
        large_ind = find(amp_diff < min_tuning/2);%1/2 max
    else
        [max_tuning,max_tuning_ind] = max(amp_diff(start_window:end_window));
        max_tuning_ind = max_tuning_ind+start_window-1;
        large_ind = find(amp_diff > max_tuning/2);%1/2 max
    end
    gaps = findgaps(large_ind);
    window = [];
    for g =1:size(gaps,1)
        gp = gaps(g,:);
        gp(gp == 0) = [];
        if any(gp == max_tuning_ind)
            window = gp;
            break
        end
    end
    if ~isempty(window)
        all_amplitude_windows{unit} = window;
    else
         all_amplitude_windows{unit} = start_window:end_window;
    end
    
    %---Modify Window if Necessary based on Fixation Durations and Remove Fixation Durations that are too Short---%
    window_width = length(all_amplitude_windows{unit});
    window_end = all_amplitude_windows{unit}(end);
    window_start = all_amplitude_windows{unit}(1);
    if window_end > end_window
        window_end = window_end-twinad1;
        median_fixation_dur = median(sac_fix_ends);
        %if more than half the fixations are shorter than window length
        %cut window down to median fixation duration othwerise we will just
        %cut fixations shorter than window end
        if median_fixation_dur < window_end
            window = all_amplitude_windows{unit};
            window(window > median_fixation_dur+twinad1) = [];
            if length(window) < min_window_width;
                all_amplitude_windows{unit} = [];
            else
                all_amplitude_windows{unit} = window;
            end
            window_end = all_amplitude_windows{unit}(end)-twinad1;
        end
        %remove fixations shorter than window as these could be
        %contaminated by the next saccade
        fixations_too_short = find(sac_fix_ends < window_end);
        sac_amps(fixations_too_short) = [];
        sac_aligned_amp(fixations_too_short,:) = [];
        sac_fix_ends(fixations_too_short) = [];
        sac_fix_starts(fixations_too_short) = [];
    end
    if window_start < twinad1-min_fix_dur
        fixations_too_short = find((twinad1+sac_fix_starts) > window_start);
        sac_amps(fixations_too_short) = [];
        sac_aligned_amp(fixations_too_short,:) = [];
        sac_fix_ends(fixations_too_short) = [];
        sac_fix_starts(fixations_too_short) = [];
    end

    %---Saccade Amplitude Analysis---%
    spike_counts = sum(sac_aligned_amp(:,all_amplitude_windows{unit}),2);
    if sum(spike_counts > 0) >  min_saccades_with_spikes*size(sac_aligned_amp,1)
        firing_rates = sum(sac_aligned_amp(:,all_amplitude_windows{unit}),2)*1000/window_width;
        [amps,amplitude_binned_firing_rate_curves{unit}] =  amplitude_bin_firing(sac_amps,firing_rates,bin_amplitude,max_amplitude);
        shuffled_corrs = NaN(1,numshuffs); %shuffled correlations
        parfor shuff = 1:numshuffs
            ind = randperm(length(sac_amps)); %shuffle indeces
            shuffled_sac_amps = sac_amps(ind); %shuffled saccade amplitudes
            [amps2,bfrc] = amplitude_bin_firing(shuffled_sac_amps,firing_rates,bin_amplitude,max_amplitude);
            shuffled_corrs(shuff) = corr(bfrc(2:end)',amps2(2:end)','row','pairwise','type','Spearman');
        end
        
        amplitude_correlations(unit) = corr(amplitude_binned_firing_rate_curves{unit}',amps','row','pairwise','type','Spearman');%observed value
        shuffled_amplitude_correlations{unit} = shuffled_corrs; %shuffled correlations
        amplitude_correlations_percentile(unit)= 100*sum(abs(amplitude_correlations(unit)) >  abs(shuffled_amplitude_correlations{unit}))/numshuffs;
    else
        firing_rates = sum(sac_aligned_amp(:,all_amplitude_windows{unit}),2)*1000/window_width;
        [amps,amplitude_binned_firing_rate_curves{unit}] =  amplitude_bin_firing(sac_amps,firing_rates,bin_amplitude,max_amplitude);
        amplitude_correlations(unit) = corr(amplitude_binned_firing_rate_curves{unit}',amps','row','pairwise','type','Spearman');%observed value
        disp('Too sparse')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Saccade Amplitude Plots---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %----Saccade Amplitude---%    
    figure
    
    %---Saccade Aligned Firing Rate Curve--%
    subplot(2,3,1)
    dofill(tm,sac_aligned_amp,'black',1,smval); %smoothed fiirng rate curve
    hold on
    yl = ylim;
    if yl(1) < 0
        ylim([0 yl(2)]);
        yl(1) = 0;
    end
    if ~isempty(all_amplitude_windows{unit})
        h = fill([all_amplitude_windows{unit}(1) all_amplitude_windows{unit}(end) ...
            all_amplitude_windows{unit}(end) all_amplitude_windows{unit}(1) all_amplitude_windows{unit}(1)]-twinad1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
        set(h,'facealpha',.25,'EdgeColor','None')
    end
    hold off
    xlim([-twinad1 twinad2]);
    xlabel('Time from Saccade Start (ms)')
    ylabel('Firing Rate (Hz)')
    title('All Saccade Aligned Activity')
    yls(1,:) = ylim;
    
    %---Plot Firing Rate Curves for Large vs Small Amplitude Saccades---%
    small = prctile(sac_amps,20);
    medium1 = prctile(sac_amps,40);
    medium2 = prctile(sac_amps,60);
    large = prctile(sac_amps,80);
    subplot(2,3,2)
    hold on
    dofill(tm,sac_aligned_amp(sac_amps < small,:),'black',1,smval); %smoothed fiirng rate curve
    dofill(tm,sac_aligned_amp(sac_amps < medium2 & sac_amps > medium1,:),'blue',1,smval); %smoothed fiirng rate curve=
    dofill(tm,sac_aligned_amp(sac_amps > large,:),'green',1,smval); %smoothed fiirng rate curve
    yl = ylim;
    if yl(1) < 0
        ylim([0 yl(2)]);
    end
    hold off
    xlim([-twinad1 twinad2]);
    xlabel('Time from Saccade Start (ms)')
    ylabel('Firing Rate (Hz)')
    legend('Small','Medium','Large')
    title(['Smallest/Medium/Largest 20% Saccades']);
    
    %---Plot Firing Rate Curves for Saccades of different amplitudes---%
    subplot(2,3,3)
    plot(tm,crude_amplitude_amplitude_curves')
    xlim([-twinad1 twinad2]);
    xlabel('Time From Saccade Start (ms)')
    ylabel('Firing Rate')
    title(['Firing Rate Curves by Sac. Ampl. in 4 dva bins']);
    box off
    
    %---Fixation Aligned Raster Sorted by Saccade amplitude for out2out & in2in fixations---%
    subplot(2,3,4)
    [~,si] = sort(sac_amps);
    sac_aligned_sorted = sac_aligned_amp(si,:);
    [trial,time] = find(sac_aligned_sorted == 1);
    plot(time-twinad1,trial,'.k')
    if ~isempty(trial)
        yl = [0 max(trial)+1];
        ylim([0 max(trial)+1]);
        ylim(yl);
        hold on
        if ~isempty(all_amplitude_windows{unit})
            h = fill([all_amplitude_windows{unit}(1) all_amplitude_windows{unit}(end) ...
                all_amplitude_windows{unit}(end) all_amplitude_windows{unit}(1) all_amplitude_windows{unit}(1)]-twinad1,...
                [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
            set(h,'facealpha',.25,'EdgeColor','None')
        end
        hold off
    end
    xlim([-twinad1 twinad2]);
    ylabel('Ranked Sac. Amp.')
    xlabel('Time from Saccade Start (ms)')
    title('Saccade Aligned Raster-Saccade amplitude')
    box off
    
    %---Firing Rate Curve for Saccade Amplitude---%
    subplot(2,3,5)
    plot(amps,amplitude_binned_firing_rate_curves{unit})
    xlabel('Saccade Ampltiude (dva)')
    ylabel('Firing Rate (Hz)')
    title(sprintf(['\\rho_{amp} = '  num2str(amplitude_correlations(unit),2) ' (' num2str(amplitude_correlations_percentile(unit),3) '%%)']))
    box off
    
    %---Firing Rate Map---%
    subplot(2,3,6)
    firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
    maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
    h = imagesc(firing_rate_map);
    set(h,'alphadata',~isnan(firing_rate_map));
    axis off
    axis equal
    colormap('jet')
    colorbar
    clim = caxis;
    caxis([clim(1) maxfr])
    title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(unit) > 95;
        title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
            '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
    end
    if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
        title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
    end
    if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
        title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
            '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
    end
    title(sprintf(title_str));
    
    
    if multiunit(unit)
        subtitle(['Multiunit ' task_file(1:8) ' ' unit_stats{1,unit}]);
    else
        subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
    end
    %%
    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_Saccade_amplitude_Analysis']);
    
    %%
    
end


save([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'],...
    'min_fix_dur','min_sac_amp','numshuffs','bin_deg','bin_deg2',...
    'mrls','uniformity_pvalue','unit_stats','smval','bin_amplitude','bin_amplitude2',...
    'max_amplitude','saccade_aligned_firing','saccade_directions','saccade_amplitudes',...
    'sac_in_out','uniformity_pvalue','direction_binned_firing_rate_curves','all_direction_windows',...
    'amplitude_correlations','amplitude_correlations_percentile',...
    'amplitude_binned_firing_rate_curves','all_amplitude_windows','fixation_starts','fixation_ends');

end

function [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,saccade_aligned_firing,time_window,saccade_directions)
%bin data into small degree bins
degrees = [0:bin_deg:360]-180;
binned_firing_rate = cell(1,length(degrees));
for bins = 2:length(degrees)
    these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
    binned_firing_rate{bins} = 1000*nansum(saccade_aligned_firing(these_dirs,time_window),2)/length(time_window);
end
degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];
mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
mrl = circ_r(saccade_directions*pi/180,sum(saccade_aligned_firing(:,time_window),2)); %MRL for unbinned data

%NOT rechecked for bugs since currently not used
%binned data
% degrees = [0:3:360]-180; %binned degrees
% binned_firing_rate = cell(1,length(degrees));
% for bins = 2:length(degrees)
%     these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
%     binned_firing_rate{bins} = 1000*nansum(saccade_aligned_firing(these_dirs,time_window),2)/length(time_window);
% end
% degrees = degrees(2:end);
% degrees = degrees*pi/180;
% degrees = [degrees degrees(1)];
% mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
% mbfr = mean_binned_firing_rate;
% mbfr(isnan(mbfr)) = 0;
%mrl = circ_r(degrees(1:end-1),mbfr,bin_deg/180*pi); %calculate for binned data
end

function [prefered_dirs,anti_prefered_dirs,smoothed_firing] = select_prefred_indeces(binned_firing,binned_directions,dirs,smval,bin_size)
orginal_dirs = dirs; %store for later since may be modifying dirs variable
bf = [binned_firing(end-(3*smval):end) binned_firing binned_firing(1:3*smval)];%buffer so don't get edege artifacts from filtering
binned_firing = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
binned_firing = binned_firing(3*smval+2:end-3*smval);%remove buffers
smoothed_firing = binned_firing;%for output

%---Calculate Prefered Direction---%
%over smoothed so small peaks don't interfer detecting prefered direction
binned_firing2 = nandens(bf,3*smval,'gauss',1);
binned_firing2 = binned_firing2(3*smval+2:end-3*smval);%remove buffers
prefered_index = find((binned_firing2(1:end-1) == max(binned_firing2(1:end-1))));
prefered_index = prefered_index(1);
prefered_direction = 180/pi*binned_directions(prefered_index);%highest firing rate

%---Calculate Anti-prefered Direction---%
%over smoothed so small peaks don't interfer detecting prefered direction
binned_firing2 = nandens(bf,3*smval,'gauss',1);
binned_firing2 = binned_firing2(3*smval+2:end-3*smval);%remove buffers
anti_index = find(binned_firing2(1:end-1) == min(binned_firing2(1:end-1)));
anti_index = anti_index(1);
anti_prefered_direction = 180/pi*binned_directions(anti_index); %lowest firing rate

%---Calculate FWHM for Prefered Direction---%
binned_firing2 = binned_firing-mean(binned_firing);
if abs(prefered_direction) > 90
    %rotate axis if prefered direction is near window edge
    zeroind = find(binned_directions == 0);
    binned_firing2 = [binned_firing2(zeroind:end) binned_firing2(1:zeroind-1)];
    
    binned_firing3 = nandens(bf,3*smval,'gauss',1);
    binned_firing3 = binned_firing3(3*smval+2:end-3*smval);%remove buffers
    binned_firing3 = [binned_firing3(zeroind:end) binned_firing3(1:zeroind-1)];
    prefered_index = find((binned_firing3(1:end-1) == max(binned_firing3(1:end-1))));
    prefered_index = prefered_index(1);
end

binned_firing2 = binned_firing2-mean(binned_firing2);
large_ind = find(binned_firing2 > binned_firing2(prefered_index)/2);%1/2 the max
gaps = findgaps(large_ind);
window = [];
for g =1:size(gaps,1)
    gp = gaps(g,:);
    gp(gp == 0) = [];
    if any(gp == prefered_index)
        window = gp;
        break
    end
end
deg_win = length(window)/2*bin_size;
if deg_win < 15
    deg_win = 15;
elseif deg_win > 45
    deg_win = 45;
end

%---Calculate FWHM for Anti-Prefered Direction---%
if abs(anti_prefered_direction) > 90
    %rotate axis if anti-prefered direction is near window edge
    zeroind = find(binned_directions == 0);
    binned_firing2 = [binned_firing2(zeroind:end) binned_firing2(1:zeroind-1)];
    
    binned_firing3 = nandens(bf,3*smval,'gauss',1);
    binned_firing3 = binned_firing3(3*smval+2:end-3*smval);%remove buffers
    binned_firing3 = [binned_firing3(zeroind:end) binned_firing3(1:zeroind-1)];
    anti_index = find(binned_firing3(1:end-1) == min(binned_firing3(1:end-1)));
    anti_index  = anti_index(1);
end

binned_firing2 = binned_firing2-mean(binned_firing2);
small_ind = find(binned_firing2 < binned_firing2(anti_index)/2);%1/2 the min
gaps = findgaps(small_ind);
window = [];
for g =1:size(gaps,1)
    gp = gaps(g,:);
    gp(gp == 0) = [];
    if any(gp == anti_index)
        window = gp;
        break
    end
end
deg_win_anti = length(window)/2*bin_size;
if deg_win_anti < 15
    deg_win_anti = 15;
elseif deg_win_anti > 45
    deg_win_anti = 45;
end

%---Select Directions in Prefered Direction---%
if prefered_direction > 180-deg_win %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif prefered_direction < -180+deg_win
    prefered_direction = prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
prefered_dirs = find((dirs <= prefered_direction+deg_win)& (dirs >= prefered_direction-deg_win));

%---Select Directions in Anti-Prefered Direction---%
dirs = orginal_dirs;%restore dirs variable since may have just modifid above
if anti_prefered_direction > 180-deg_win_anti %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif anti_prefered_direction < -180+deg_win_anti
    anti_prefered_direction = anti_prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
anti_prefered_dirs = find((dirs <= anti_prefered_direction+deg_win_anti)& (dirs >= anti_prefered_direction-deg_win_anti));
end


function [amps,binned_firing_rates] = amplitude_bin_firing(sacamps,firing_rates,bin_amplitude,max_amplitude)
amps = [2:bin_amplitude:max_amplitude];
binned_firing_rates = NaN(1,length(amps));
for bin = 1:length(amps)-1
    these_amplitudes = (sacamps < amps(bin+1) & sacamps >= amps(bin));
    binned_firing_rates(bin)= mean(firing_rates(these_amplitudes));
end
end

function [max_tuning_ind,window] = get_half_max_window(values,tmin,tmax)
if ~isempty(tmin) && ~isempty(tmax)
    [max_tuning,max_tuning_ind] = max(values(tmin:tmax));
    max_tuning_ind = max_tuning_ind+tmin-1;
else
    [max_tuning,max_tuning_ind] = max(values);
end
large_ind = find(values > max_tuning/2);%1/2 the max
gaps = findgaps(large_ind);
window = [];
for g =1:size(gaps,1)
    gp = gaps(g,:);
    gp(gp == 0) = [];
    if any(gp == max_tuning_ind)
        window = gp;
        break
    end
end
end
