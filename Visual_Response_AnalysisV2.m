function Visual_Response_AnalysisV2(data_dir,figure_dir,session_data)
% written by Seth Konig July 11, 2016 updated to V2 August 31, 2016 to focus on
% visual response and time delay so can look at lag for eye movements and
% events. Code determines if neurons show consistent, peaky response to
% cross hair (& fixation on cross hair), image on for the first 1 second,
% and image on for the first 5 seconds of response, and image off.
%
% code rechecked for bugs on January 15, 2017 by SKD

figure_dir = [figure_dir 'Visual Response\'];
twin1 = 200;% how much time to take before event cross appears and how much to ignore during fixation
twin2 = 1000;%how much time to look at after stimulus onset for short window
twin3 = 500;%minimum fixation on image cross hair duraiton
twin4 = 5000; %for long window on image on
numshuffs = 10000; %number of shuffles to do for bootstrapping
smval =60; %30 ms std, 1/2 width of gaussian smoothing filter
smval2 = 300;%150 ms std, 1/2 width of gaussian smoothing filter

% time_locked_firing since ITI period is defined by 2 events (15 & 16)
event_codes = [35 8 23 24 23];
trial_start_code = 15;
task = 'ListSQ';
min_blks = 2;

sliding_window = twin1; %window width of ~200 ms
sliding_step = 25;  %step size in ms, want to be relatively small compared to sliding window wiidth
p_thresh = 0.01;%significance level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0
    return %if no units exit function
end
num_trials = length(cfg.trl); %number of trials

%get trials in which spikes are valid i.e. trials for which neuron is
%stable for
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%get important task specific information
[itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import & reformat data so that spikes are locked to events---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Aligning spike times to trial events')

%preallocate space and parallel structure of cfg
time_lock_firing = cell(num_units,length(event_codes));%event aligned spike trains
which_images = cell(1,num_units); %image #
trial_number = cell(1,num_units); %# in cfg so can import fixstats at later time
nvr = cell(1,num_units);%novel vs repeat
fix_rt = cell(1,num_units); %how long does it take on average to fixate cross once displayed
for unit = 1:num_units
    time_lock_firing{unit,1} = NaN(192,twin1+twin3);%cross hair on
    time_lock_firing{unit,2} = NaN(192,twin1+twin2);%fix crosshair
    time_lock_firing{unit,3} = NaN(192,twin1+twin2);%imgon
    time_lock_firing{unit,4} = NaN(192,2*twin3);%imgoff.
    time_lock_firing{unit,5} = NaN(192,twin1+twin4);%img on long window
    fix_rt{unit}=  NaN(1,192); %how long does it take on average to fixate cross once displayed
    nvr{unit} = NaN(1,length(cfg.trl)); %novel vs repeat
    which_images{unit} = NaN(1,length(cfg.trl)); %image #
    trial_number{unit} = NaN(1,length(cfg.trl)); %# in cfg so can import fixstats at later time
end

for t = 1:num_trials
    if any(cfg.trl(t).allval == event_codes(3)); %in which image was displayed or 1st item in sequence was displayed
        if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(end) %then sequence trial
            continue % go to the next trial
        end
        
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
        trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after image off
        crosson = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(1))-trial_start; %cross hair on
        fix_cross = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(2))-trial_start; %fixation on cross hair according to cortex
        fix_cross = fix_cross(1);%since fixation on image also counts as event 8
        imgon = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(3))-trial_start; %when image turns on
        imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(4))-trial_start; %when image turns off
        
        if imgon-fix_cross > 750%cortex sometimes has minimal lag up to ~20 ms
            imgon2 =  fix_cross+750;
        else
            imgon2 = imgon;
        end
        
        for unit = 1:num_units
            if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                
                %---image info---%
                img_index = find(cfg.trl(t).cnd == img_cnd); %image index
                if any(isnan(which_img(img_index)))
                    continue
                end
                nvr{unit}(img_index) = novel_vs_repeat(img_index); %whether image was novel or repeated
                which_images{unit}(img_index) = which_img(img_index); %image number
                fix_rt{unit}(img_index) = fix_cross-crosson; %reaction time to fixate on cross
                trial_number{unit}(img_index) = t; %trial #
                
                spikes = find(data(unit).values{t}); %spike trains for this trial
                
                %---crosson locked firing---%
                event_spikes = spikes(spikes > crosson-twin1 & spikes <= crosson+twin3)-crosson+twin1;
                tempvec = zeros(1,twin1+twin3);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,1}(img_index,:) = tempvec;
                
                %---fix cross locked firing---%
                %variable duration event
                event_spikes = spikes(spikes > fix_cross-twin1 & spikes <= imgon2)-fix_cross+twin1;
                timevec = [fix_cross-twin1+1 imgon2]-fix_cross+twin1;
                tempvec = zeros(1,timevec(2)-timevec(1)+1);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,2}(img_index,timevec(1):timevec(2)) = tempvec;
                
                %---image on locked firing---%
                event_spikes = spikes(spikes > imgon-twin1 & spikes <= imgon+twin2)-imgon+twin1;
                tempvec = zeros(1,twin1+twin2);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,3}(img_index,:) = tempvec;
                
                %---image off locked firing---%
                if t ~= length(cfg.trl)
                    %will need to look in next trial to get rest of data
                    %since usually only get ~200 ms of data
                    rest_time = twin3-(trial_end-imgoff);%twin3-time from image off to trial end
                    spikes2 = find(data(unit).values{t+1}); %next trial's spike trains
                    spikes2(spikes2 > rest_time) = [];%remove next trial's spike trains spikes up to rest of time
                    spikes2 = spikes2+(twin3-rest_time)+twin3;%reindex
                    event_spikes = spikes(spikes > imgoff-twin3 & spikes <= imgoff+twin3)-imgoff+twin3;%this trials spikes
                    event_spikes = [event_spikes spikes2]; %concatenate this plus next trial's spikes
                    tempvec = zeros(1,2*twin3);
                    tempvec(event_spikes) = 1;
                    time_lock_firing{unit,4}(img_index,:) = tempvec;
                end
                
                %---long image on locked firing---%
                event_spikes = spikes(spikes > imgon-twin1 & spikes <= imgon+twin4)-imgon+twin1;
                tempvec = zeros(1,twin1+twin4);
                tempvec(event_spikes) = 1;
                time_lock_firing{unit,5}(img_index,:) = tempvec;
                
            end
        end
    end
end

%remove excess NaNs associated with error trials
time_lock_firing = laundry(time_lock_firing);
nvr = laundry(nvr);
which_images = laundry(which_images);
fix_rt = cellfun(@nanmean,fix_rt);
trial_number = laundry(trial_number);

%for comparing novel to repeat only want the trials in which both
%images were shown
for unit = 1:num_units
    rmv = []; %images to remove
    for img = 1:96
        ind = find(which_images{unit} == img);
        if length(ind) ~= 2 %so either novel or repeat but not both
            rmv = [rmv ind];
        end
    end
    which_images{unit}(rmv) = NaN; %set to NaN to keep structure
    nvr{unit}(rmv) = NaN; %set to NaN to keep structure
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Perform Temporal analysis and signifance testing---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Determining if spikes are locked to trial events')
Fs = data(1).fsample; %should be 1000
info_type = 'temporal';

epoch_data = [];
epoch_data.rate = NaN(num_units,size(time_lock_firing,2));
epoch_data.temporalstability = NaN(num_units,size(time_lock_firing,2));
epoch_data.shuffled_rate = cell(num_units,size(time_lock_firing,2));
epoch_data.shuffled_temporalstability = cell(num_units,size(time_lock_firing,2));
epoch_data.rate_prctile = NaN(num_units,size(time_lock_firing,2));
epoch_data.temporalstability_prctile = NaN(num_units,size(time_lock_firing,2));
for unit = 1:num_units
    for event = 1:length(event_codes);
        if isempty(time_lock_firing{unit,event})
            continue
        elseif nansum(nansum(time_lock_firing{unit,event})) == 0; %no spikes
            epoch_data.rate(unit,event) = 0;
            epoch_data.rate_prctile(unit,event) = 0;
            epoch_data.temporalstability_prctile(unit,event) = 0;
        else
            if event == 5 %long image period only difference is smval2
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval2,Fs);
            elseif event == 2 %fixation on cross variable duration until image is on so only take first part up to twin3
                if nansum(nansum(time_lock_firing{unit,event}(:,1:twin3+twin1))) > 0
                    [temporal_info,shuffled_info_rate] = ...
                        estimated_mutual_information(time_lock_firing{unit,event}(:,1:twin1+twin3),numshuffs,info_type,smval,Fs);
                else
                    epoch_data.rate(unit,event) = 0;
                    epoch_data.rate_prctile(unit,event) = 0;
                    epoch_data.temporalstability_prctile(unit,event) = 0;
                end
            else
                [temporal_info,shuffled_info_rate] = ...
                    estimated_mutual_information(time_lock_firing{unit,event},numshuffs,info_type,smval,Fs);
            end
            
            epoch_data.rate(unit,event) = temporal_info.skaggs;
            %only want 1st vs 2nd half no e/o since not using
            epoch_data.temporalstability(unit,event) = temporal_info.temporalstability(1);
            epoch_data.shuffled_rate{unit,event} = shuffled_info_rate.skaggs;
            epoch_data.shuffled_temporalstability{unit,event} = shuffled_info_rate.temporalstability(1,:);
            epoch_data.rate_prctile(unit,event) = 100*sum(temporal_info.skaggs > shuffled_info_rate.skaggs)/numshuffs;
            epoch_data.temporalstability_prctile(unit,event) = 100*sum(temporal_info.temporalstability(1) >...
                shuffled_info_rate.temporalstability(1,:))/numshuffs;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Determine if Neuron is Visually Responsive using Method Similar to Mike Jutras---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%i.e. sliding window analysis
sig_visual_response = zeros(num_units,twin1+twin4); %for short and long window
baseline_firing_rate = NaN(2,num_units); %mean and std
p_visual_response = NaN(num_units,201); %p value by step
for unit = 1:num_units
    if ~isempty(time_lock_firing{unit,5})%same as time_lock_firing{unit,3) just full window
        
        
        %define baseline firing rate of neuron as 200 ms before image on
        %in case there is a fixation response or crosshair response
        
        %twin1 get to event start then want to go another twin1 to ignore 200 ms
        avg_firing_rate = 1000*nansum(time_lock_firing{unit,5}(:,1:twin1),2)./twin1; %trial-by-trial firing rate
        baseline_firing_rate(1,unit) = mean(avg_firing_rate);
        baseline_firing_rate(2,unit) = std(avg_firing_rate);
        fr_distribution = sum(time_lock_firing{unit,5}(:,1:twin1),2); %trial-by-trial spike count distribution
        
        %---Analysis of Visual Response for Short Window---%
        total_time = size(time_lock_firing{unit,5},2);
        for step = 1:(total_time/sliding_step)-sliding_window/sliding_step+1
            ind = sliding_step*(step-1)+1:sliding_window+sliding_step*(step-1);
            [~,p_visual_response(unit,step)] = kstest2(sum(time_lock_firing{unit,5}(:,ind)'),...
                fr_distribution);
            if p_visual_response(unit,step) < p_thresh
                sig_visual_response(unit,ind) = 2;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t12 = -twin1:twin2-1;
t13 = -twin1:twin3-1;
t133 = -twin1:twin3-1+250;
t3 = -twin3:twin3-1;
t14 = -twin1:twin4-1;
unit_names = unit_stats(1,:);

for unit = 1:num_units
    if ~isempty(time_lock_firing{unit,1})
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Peri Image Event Plots---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        
        
        %---Plot Firing rate curve locked to Cross and Fixation on Cross---%
        subplot(2,2,1)
        hold on
        [~,~,~,y] =dofill(t13-round(fix_rt(unit)),time_lock_firing{unit,1},'black',1,smval);%smoothed cross onset curve
        enough_samples = find(sum(~isnan(time_lock_firing{unit,2})) > 40); %variable event length so cut time points without enough samples
        [~,~,~, y1,~] =dofill2(enough_samples-twin1,time_lock_firing{unit,2}(:,enough_samples),'green',1,smval); %smoothed fixation on cross curve
        plot([-twin3 twin3+250],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--') %pre-image "baseline" line
        hold off
        xlim([-twin3 twin3+250])
        set(gca,'Xtick',[-500 -250 0 250 500 750])
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Fixation on Cross (ms)')
        legend('Cross On','Fix on Cross','"Baseline"')
        if epoch_data.rate_prctile(unit,2) > 90 || epoch_data.temporalstability_prctile(unit,2) > 90
            title(['bit_{fix} ' num2str(epoch_data.rate_prctile(unit,2),3) '% ' ...
                '\rho_{fix} = ' num2str(epoch_data.temporalstability(unit,2),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,2),3)  '%' ])
        end
        box off
        
        %---Plot Raster aligned to fixation on cross---%
        subplot(2,2,3)
        [trial,time] = find(time_lock_firing{unit,2} == 1);
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin3 twin3+250])
            ylim([0 max(trial)+1]);
        end
        set(gca,'Xtick',[-500 -250 0 250 500 750])
        ylabel('Trial #')
        xlabel('Time from Fixation on Cross (ms)')
        box off
        
        
        %---Plot Firing rate curve locked to Image Off---%
        subplot(2,2,2)
        [~,~,~,y2] =dofill(t3,time_lock_firing{unit,4},'black',1,smval); %smoothed image off firing rate curve
        xlim([-twin3 twin3])
        set(gca,'Xtick',[-500 -250 0 250 500])
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image Off (ms)')
        if epoch_data.rate_prctile(unit,4) > 90 || epoch_data.temporalstability_prctile(unit,4) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,4),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,4),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,4),3)  '%' ])
        end
        box off
        
        %---Plot Raster aligned to Image off---%
        subplot(2,2,4)
        [trial,time] = find(time_lock_firing{unit,4} == 1);
        if ~isempty(trial)
            plot(time-twin3,trial,'.k')
            xlim([-twin3 twin3])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image Off (ms)')
        set(gca,'Xtick',[-500 -250 0 250 500])
        box off
        
        
        %---Scale plots to be the same y-axis---%
        ymin = 0.85*min([y y1 y2]);
        if ymin < 0.1
            ymin = 0;
        end
        ymax = 1.2*max([y y1 y2]);
        if ymin ~= ymax
            subplot(2,2,1)
            ylim([ymin ymax])
            subplot(2,2,2)
            ylim([ymin ymax])
        end
        
        n_str = [' Fixation on Cross and Image Off, n_ =' num2str(size(time_lock_firing{unit,1},1))];
        if multiunit(unit)
            subtitle(['Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['' unit_names{unit} n_str]);
        end
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Visual Response_CrossHair']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Image Onset Event Plots---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure
        
        %---Plot Firing rate locked to Image On 1 second window---%
        subplot(2,2,1)
        hold on
        [~,~,~,y3] =dofill(t12,time_lock_firing{unit,3},'black',1,smval); %smoothed image onset curve
        plot([-twin1 twin2],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')%pre-image "baseline" line
        hold off
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image On (ms)')
        xlim([-twin1 twin2])
        if epoch_data.rate_prctile(unit,3) > 90 || epoch_data.temporalstability_prctile(unit,3) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
        end
        box off
        
        %---Plot raster for image on 1 second window---%
        %in case we want to look at latency or something
        subplot(2,2,3)
        [trial,time] = find(time_lock_firing{unit,3} == 1); %novel
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin2])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Firing rate locked to Image On 5 second window---%
        subplot(2,2,2)
        [~,~,~,y4] =dofill(t14,time_lock_firing{unit,5},'black',1,smval2); %smoothed firing rate curve
        hold on
        plot([-twin1 twin4],[baseline_firing_rate(1,unit) baseline_firing_rate(1,unit)],'k--')%pre-image "baseline" line
        hold off
        ylabel('Firing Rate (Hz)')
        xlabel('Time from Image On (ms)')
        xlim([-twin1 twin4])
        if epoch_data.rate_prctile(unit,5) > 90 || epoch_data.temporalstability_prctile(unit,5) > 90
            title(['bit ' num2str(epoch_data.rate_prctile(unit,5),3) '% ' ...
                '\rho = ' num2str(epoch_data.temporalstability(unit,5),2) ' ' ...
                num2str(epoch_data.temporalstability_prctile(unit,5),3)  '%' ])
        end
        box off
        
        %---Raster for Image On 5 second window---%
        subplot(2,2,4)
        [trial,time] = find(time_lock_firing{unit,5} == 1);
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin4])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Scale plots to be the same and plot sig. difference from pre-image "baseline"---%
        ymin = 0.85*min([y3 y4]);
        if ymin < 0.1
            ymin = 0;
        end
        ymax = 1.2*max([y3 y4]);
        if ymin ~= ymax
            
            %short 1 second window
            subplot(2,2,1)
            hold on
            ylim([ymin ymax])
            gaps = findgaps(find(( sig_visual_response(unit,1:twin1+twin2))));
            if ~isempty(gaps)
                for g = 1:size(gaps,1)
                    gp = gaps(g,:);
                    gp(gp == 0) = [];
                    h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                        [ymin ymin ymax ymax ymin],'k');
                    uistack(h,'down')
                    set(h,'facealpha',.25,'EdgeColor','None')
                end
            end
            hold off
            
            %long 5 second window
            subplot(2,2,2)
            hold on
            ylim([ymin ymax])
            gaps = findgaps(find((sig_visual_response(unit,:))));
            if ~isempty(gaps)
                for g = 1:size(gaps,1)
                    gp = gaps(g,:);
                    gp(gp == 0) = [];
                    h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                        [ymin ymin ymax ymax ymin],'k');
                    uistack(h,'down')
                    set(h,'facealpha',.25,'EdgeColor','None')
                end
            end
            hold off
        end
        
        n_str = [' Image Onset, n_ =' num2str(size(time_lock_firing{unit,1},1))];
        if multiunit(unit)
            subtitle(['Multiunit ' unit_names{unit} n_str]);
        else
            subtitle(['' unit_names{unit} n_str]);
        end
        
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Image_Visual Response']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'time_lock_firing',...
    'smval','smval2','epoch_data','unit_stats','twin1','twin2','twin3','twin4',...
    'p_thresh','sliding_window','sliding_step','baseline_firing_rate',...
    'sig_visual_response','p_visual_response','unit_names','nvr','which_images',...
    'fix_rt','trial_number')
disp(['Time Locked Data Analyis for ' task_file(1:8) ' saved']);
end