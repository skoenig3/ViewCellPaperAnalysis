function List_Fixation_Analysis(data_dir,figure_dir,session_data)
% Code only looks at fixation aligned data instead of saccade aligned.
%
% Function analyizes spike times correlated with eye movements in the
% List image portion of the ListSQ task.
%
% Inputs:
%   1) data_dir: directory where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: contains relavent session information
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-Eyemovement_Locked_List_results.mat'
%
% Code rechecked for bugs on January 10-11,2017 SDK

figure_dir = [figure_dir 'List Fixation Analysis\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important Analysis parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin1 = 200;% how much time to take before eye movement starts, in case neurons are encoding upcomming eye movement
twin2 = 400;%how much time to take after eye movement has started
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer
trial_start_code = 15; %start of ITI/trial
trial_end_code = 20;
imgon_code = 23; %image turned on
imgoff_code = 24; %image turned off
smval = 30 ;%15 ms std, want to smooth at high frequency modulations
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to determine significance
info_type = 'temporal';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)

min_fix_dur = 100; %100 ms, don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva, don't want mini/micro saccades too small and hard to detect
min_num_fix = 250; %at least 250 fixatoins with a certain duration to analyze for time limited to fixation duration

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

%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end
Fs = data(1).fsample; %should be 1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Get successful trials---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([task_file ': running Fixation modulation anlaysis..'])

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,imgon_code);


%preallocate space and parallel structure of cfg
num_trials = length(cfg.trl); %total number of trials
image_trials = zeros(1,num_trials); %image trials only
for t = 1:num_trials %only take trials in which image was actually shown
    if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
        image_trials(t) = 1;
    end
end

%remove excess data associated with non image trials
image_trials = find(image_trials);
num_trials = length(image_trials); %total number of image trials
fixationstats = fixationstats(image_trials); %only take the eye data from successful trials


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Process eye data locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixation_info = NaN(9,25*length(fixationstats));
%1) trial #
%2) %fixation start time from trial start
%3) ordinal fixation #
%4) fixation start relative to image onset
%5) fixation duration
%6) preceeding saccade amplitude
%7) preceeding saccade duration
%8) preceeding saccade direction 
%9) image novel/repeat

fix_ind = 1; %start fixation index @ 1
for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
    fixationtimes = fixationstats{t}.fixationtimes; %start and end times of fixations 
    saccadetimes = fixationstats{t}.saccadetimes; % start and end times of saccades
    xy = fixationstats{t}.XY; %x and y eye data 
    
    trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code); %time at start of trial 
    trial_end = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_end_code);%time at end of trial
    imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgon_code)-trial_start; %time of image on relative to start of trial
    imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == imgoff_code)-trial_start; %time of image off relative to start of trial
    nvr = novel_vs_repeat(img_cnd == cfg.trl(image_trials(t)).cnd); %novel or repeat image
    
    if any(isnan(nvr))%cortex presentation error so skip trial, see get_image_numbers.m
        continue
    end
    
    % if monkey isn't paying attention data isn't probably
    % worth much plus have to cut off somewhere
    if imgoff-imgon > 1.5*imgdur-1
        imgoff = imgon+1.5*imgdur-1;
    end
    
    %find fixations and saccades that did not occur during the image period;
    %should also take care of the 1st fixation on the crosshair
    
    %fixation started before image turned on
    invalid= find(fixationtimes(1,:) < imgon);
    fixationtimes(:,invalid) = [];
    
    %fixation ended after the image turned off so firing rate could corrupted by image turning off
    invalid= find(fixationtimes(2,:) > imgoff);
    fixationtimes(:,invalid) = [];
    
    %saccade started before image turned on
    invalid= find(saccadetimes(1,:) < imgon);
    saccadetimes(:,invalid) = [];
    
    %Saccaded ended after the image turned off so firing rate could corrupted by image turning off
    invalid= find(saccadetimes(2,:) > imgoff);
    saccadetimes(:,invalid) = [];
    
    for f = 1:size(fixationtimes,2);
        prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
        if isempty(prior_sac) %trial ended or eye outside of image
            continue %try next one
        end
        sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
        fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
        if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
            fixation_info(1,fix_ind) = t; %trial #
            fixation_info(2,fix_ind) = fixationtimes(1,f); %fixation start time from trial start
            fixation_info(3,fix_ind) = f; %ordinal fixation #
            fixation_info(4,fix_ind) = fixationtimes(1,f)-imgon; %fixation start relative to image onset
            fixation_info(5,fix_ind) = fix_dur;%duration of this fixation
            fixation_info(6,fix_ind) = sacamp; %preceeding saccade's amplitude
            fixation_info(7,fix_ind) = saccadetimes(2,prior_sac)-saccadetimes(1,prior_sac)+1;%preceediing saccade's duration to know how far to go back
            fixation_info(8,fix_ind) = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%saccade direction
            fixation_info(9,fix_ind) = nvr; %novel vs repeat images
            fix_ind = fix_ind+1; %add index +1 
        end
    end
end

%Remove excess NaNs;
fixation_info = laundry(fixation_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Firing Rate Locked to Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate space for variables
fixation_locked_firing = cell(1,num_units); %firing rate locked to fixation
fixation_information = cell(1,num_units); %start time relative to image onset and ordinal #
for unit = 1:num_units
    fixation_locked_firing{unit}= NaN(size(fixation_info,2),twin1+twin2);
    fixation_information{unit} = NaN(size(fixation_info,2),9); %going to organize by unit so parallels with spike structure
end

fix_ind = ones(1,num_units); %fixation index # by unit since may have different trial counts
for trial = 1:num_trials
    for unit = 1:num_units;
        if image_trials(trial) >= valid_trials(1,unit) && image_trials(trial) <= valid_trials(2,unit) %only valid trials
            fixations = find(fixation_info(1,:) == trial); %fixations for this trial 
            spikes = find(data(unit).values{image_trials(trial)}); %spikes for this trail
            
            %collect spike times relative to fixation start
            for f = 1:length(fixations);
                fixt = fixation_info(2,fixations(f)); %fixation start relative to trial start
                fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1; %spikes -twin1:twin2 around fixation
                temp = zeros(1,twin1+twin2); %zeros to put spikes in
                temp(fix_spikes) = 1;%1 when spike occured
                fixation_locked_firing{unit}(fix_ind(unit),:) = temp;
                fixation_information{unit}(fix_ind(unit),:) = fixation_info(:,fixations(f))';
                fix_ind(unit) = fix_ind(unit)+1; %+1 to index for that unit
            end
        end
    end
end
%remove excess NaNs
fixation_locked_firing = laundry(fixation_locked_firing);
fixation_information = laundry(fixation_information);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Mutual Info for Eye Movements and Firing Rate---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temporal_info.fixation = [];
temporal_info.fixation.rate = NaN(1,num_units);
temporal_info.fixation.shuffled_rate = cell(1,num_units);
temporal_info.fixation.shuffled_rate_prctile = NaN(1,num_units);
temporal_info.fixation.temporalstability = NaN(2,num_units);  %row 1 by half, row 2 even/odd trials
temporal_info.fixation.shuffled_temporalstability = cell(1,num_units);
temporal_info.fixation.shuffled_temporalstability_prctile = NaN(2,num_units);%row 1 by half, row 2 even/odd trials

for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %don't want to run on trials with the first eye movements occuring with 500 (twin) ms of image onset
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:);
        
        %since many units appear sparse especially spatial ones want to
        %run on eye movemetns that actually have spikes so use these only
        %doesn't change statistics (only scale firing rates) just decreases run time
        fixation_firing(nansum(fixation_firing,2) == 0,:) = [];%remove fixations without spikes
        
        [observed_info_rate,shuffled_info_rate]...
            = estimated_mutual_information(fixation_firing,numshuffs,info_type,smval,Fs);
        temporal_info.fixation.rate(unit) = observed_info_rate.skaggs; %skaggs information score
        temporal_info.fixation.shuffled_rate{unit} = shuffled_info_rate.skaggs; %shuffled skaggs information score
        temporal_info.fixation.shuffled_rate_prctile(unit) = 100*sum(...
            observed_info_rate.skaggs > shuffled_info_rate.skaggs)/numshuffs; %observed skaggs information score shuffled percentile
        temporal_info.fixation.temporalstability(:,unit) = observed_info_rate.temporalstability;%temporal correlation: half,even/odd
        temporal_info.fixation.shuffled_temporalstability{unit} = shuffled_info_rate.temporalstability; %shuffled temporal correlation
        temporal_info.fixation.shuffled_temporalstability_prctile(1,unit) = ...
            100*sum(observed_info_rate.temporalstability(1) > ...
            shuffled_info_rate.temporalstability(1,:))/numshuffs; %temporal correlation half shuffled percentile
        temporal_info.fixation.shuffled_temporalstability_prctile(2,unit) = ...
            100*sum(observed_info_rate.temporalstability(2) > ...
            shuffled_info_rate.temporalstability(2,:))/numshuffs;%temporal correlation even/odd shuffled percentile
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Rasters by Various Variables of Interst---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = -twin1:twin2-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ~isempty(fixation_locked_firing{unit})
        
        %---For Fixations---%
         %don't want to run on trials with the first eye movements occuring
         %with 500 (twin) ms of image onset since strong visual response in
         %many neurons
        info = fixation_information{unit}((fixation_information{unit}(:,4) > image_on_twin),:);  %fixation info for unit
        fixation_firing = fixation_locked_firing{unit}(fixation_information{unit}(:,4) > image_on_twin,:); %spike times for unit
        n_trials = round(size(fixation_firing,1)/2); %number of fixations
        
        
        figure
        
        %---plot Fixation Aligned Firing Rate curve---%
        hold on
        subplot(3,3,1)
        dofill(t,fixation_firing(1:n_trials,:),'blue',1,smval);%even trials
        dofill(t,fixation_firing(n_trials+1:end,:),'red',1,smval);%odd trials
        dofill(t,fixation_firing,'black',1,smval); %all trials
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        legend('1^{st} 1/2','2^{nd} 1/2','All','Location','Best')
        ylim(yl);
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Fixation Start (ms)')
        grid on
        ylim(yl);
        title(['Bit ' num2str(temporal_info.fixation.shuffled_rate_prctile(unit),3) '% '...
            '\rho_{1/2} = ' num2str(temporal_info.fixation.temporalstability(1,unit),2) ...
            ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(1,unit),3) ...
            '%) \rho_{e/o} = ' num2str(temporal_info.fixation.temporalstability(2,unit),2) ...
            ' (' num2str(temporal_info.fixation.shuffled_temporalstability_prctile(2,unit),3) '%)'])
        
        %---plot fixation aligned raster---%
        subplot(3,3,4)
        [trial,time] = find(fixation_firing == 1);
        plot(time-twin1,(trial),'.k')
        hold on
        plot([0 0],[0 max(trial)],'r--')
        hold off
        ylim([0 max(trial)])
        ylabel('Occurence #')
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        xlabel('Time from Fixation Start (ms)')
        title('Raster whole period')
        box off        
        
        %---Plot Fixation Aligned Firing Rate curve Limited to 1 fixation and saccade---%
        info(info(:,5) > twin2,5) = twin2; %set max fixation duration to <= twin2
        info(info(:,7) >= twin1,7) = twin1-1; %set max prior saccade duration to <= twin1
        
        %want to look only looking at window around saccade/fixation
        limited_firing = NaN(size(fixation_firing));
        info(info(:,5) > twin2,5) = twin2; %if next fixation duration is > twin set to twin
        for f = 1:size(limited_firing,1)
            ind = twin1-info(f,7):twin1+info(f,5);
            limited_firing(f,ind) = fixation_firing(f,ind);
        end
        
        %get firing rate curve limited to 1 fixation and saccade
        num_not_nans = sum(~isnan(limited_firing));%will be fewer spikes
        not_enough_samples = find(num_not_nans < min_num_fix);%.5*size(limited_firing,1)); %median duration
        [~,limited_firing_rate] = nandens3(limited_firing,smval,Fs);
        limited_firing_rate(:,not_enough_samples) = NaN;
        limited_firing_rate(:,1:twin1-44) = NaN; %don't want to go back more than median saccade duration
        limited_firing(:,1:twin1-44) = NaN; %don't want to go back more than median saccade duration
        
        subplot(3,3,2)
        plot(t,nanmean(limited_firing_rate),'k');
        hold on
        yl = ylim;
        if yl(1) < 0;
            yl(1) =0;
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        ylim(yl); %scale same as whole period firing rate
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Fixation Start (ms)')
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        title('Firing Rate Curve: 1 Fixation/Saccade')
        box off 
        
        %---Plot Raster For Limited to 1 fixation and saccade---%
        fix_durs = info(:,5);
        [~,sorted_fix_durs] = sort(fix_durs); %sort by increasing fixation duration
        %make raster plot for spikes limited time period of 1 fixation around Fixation
        lmf = limited_firing(sorted_fix_durs,:);
        subplot(3,3,5)
        [trial,time] = find(lmf == 1);
        plot(time-twin1,trial,'.k')
        hold on
        plot([0 0],[0 max(trial)],'r--')
        hold off
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylim([0 max(trial)])
        set(gca,'Xtick',[-100 0 100 200 300 400])
        ylabel('Sorted by Fixation Dur')
        xlabel('Time from Fixation Start (ms)')
        title('Raster Limited to 1 Eye Movement')
        box off 
        
        
        %---plot raster sorted by spike count---%
        f1 = sum(fixation_firing,2);
        [~,fi] = sort(f1);
        [trial,time] = find(fixation_firing(fi,:) == 1);
        subplot(3,3,3)
        plot(time-twin1,(trial),'.k')
        hold on
        plot([0 0],[0 max(trial)],'r--')
        hold off
        ylim([0 max(trial)])
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylabel('Ranked Spike Count')
        xlabel('Time from Fixation Start (ms)')
        title('For all Fixations including zero spike counts')
        box off 
        
        %---plot raster over time by Fixation order within image presentation---%
        [~,fi] = sort(info(:,3));
        subplot(3,3,6)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold on
        plot([0 0],[0 max(trial)],'r--')
        hold off
        ylim([0 max(trial)])
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylabel('Ranked Ordinal Fix #')
        xlabel('Time from Fixation Start (ms)')
        title('Discrete Time within Image Period')
        box off 
        
        %         %plot raster by Fixation amplitude
        %         [~,fi] = sort(info(:,5));
        %         subplot(3,3,7)
        %         [trial,time] = find(fixation_firing(fi,:) == 1);
        %         plot(time-twin1,trial,'.k')
        %         hold off
        %         ylim([0 max(trial)])
        %         set(gca,'Xtick',[-100 0 100 200 300 400])
        %         ylabel('Ranked Fixation Amplitude')
        %         xlabel('Time from Fixation Start (ms)')
        %         title('Fixation Amplitude')
        %           box off 
        
        %---plot raster by novel vs repeat image presentations---%
        subplot(3,3,7)
        [trial,time] = find(fixation_firing(info(:,9) == 1,:) == 1);
        plot(time-twin1,trial,'.b')
        if ~isempty(trial)
            b4 = max(trial);
        else
            b4 = 0;
        end
        hold on
        [trial,time] = find(fixation_firing(info(:,9) == 2,:) == 1);
        trial= trial+b4;
        plot(time-twin1,trial,'.r')
        if ~isempty(trial)
            ylim([0 max(trial)])
            plot([0 0],[0 max(trial)],'k--')
        end
        hold off
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylabel('Occurence #')
        xlabel('Time from Fixation Start (ms)')
        title('Novel (blue) vs Repeat (red)')
        box off 
        
        %---Plot raster by Saccade direction---%
        [~,fi] = sort(info(:,8));
        subplot(3,3,8)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        plot(time-twin1,trial,'.k')
        hold on
        plot([0 0],[0 max(trial)],'r--')
        hold off
        ylim([0 max(trial)])
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylabel('Ranked Sac. Direction')
        xlabel('Time from Fixation Start (ms)')
        title('Saccade Direction')
        box off
        
        %---plot raster by fixation duration---%
        [~,fi] = sort(info(:,5));
        subplot(3,3,9)
        [trial,time] = find(fixation_firing(fi,:) == 1);
        dur = info(fi,5);
        dur(dur > twin2) = twin2;
        plot(time-twin1,trial,'.k')
        hold on
        plot(dur,1:length(fi),'r.')
        plot([0 0],[0 max(trial)],'r--')
        hold off
        ylim([0 max(trial)])
        xlim([-twin1 twin2])
        set(gca,'Xtick',[-200 0 200 400])
        ylabel('Ranked Fixation Duration')
        xlabel('Time from Fixation Start (ms)')
        title('Fixation Duration')
        box off 
        
        n_str = [' n_{fix} =' num2str(size(fixation_locked_firing{1,unit},1))];
        if multiunit(unit)
            subtitle(['Fixation-Algined Multiunit ' task_file(1:8) ' ' unit_names{unit} n_str]);
        else
            subtitle(['Fixation-Aligned ' task_file(1:8) ' '  unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} 'Eye_Locked_analysis_Fixation_Rasters']);
    end
end

save([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat'],...
    'twin1','twin2','image_on_twin','smval','fixation_locked_firing',...
    'fixation_information','temporal_info','unit_names','image_trials','numshuffs');
end