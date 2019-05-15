function spatial_analysisV2(data_dir,figure_dir,session_data,task)
% written by Seth Konig August, 2014
% Function analyizes where the monkey is attending when spikes.
% updated SDK 1/11/16 to handlde new format and partial session data for
% vaild trials only.CVTNEW section on 1/19/16
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) figure_dir: location of where to put save figures
%   3) session_data: data containing relevant session data
%   4) task: what task was this data come from with i.e. 'cvtnew','ListSQ'
%
% Outputs:
%   1) Saves figures to figure_dir
%   2) Saves processed data to data_dir tagged with '-spatial_analysis_results.mat'
%
% Code rechecked ListSQ section for bugs October 17-18, 2016 SDK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numshuffs = 1000;% # of shuffles for bootstraspping, recommend this is at least 1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials',...
    'hdr','whole_session_mean_firing_rate','excitatory_inhibitory');

%grab unit data
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

%---determine number of blocks unit was stable for
%remove units with too few trials
%these are the absolute minimum data required to do data analysis
if strcmpi(task,'cvtnew')
    min_blks = 2; %~100 trials may want to change
else
    min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
    %this means at least 64 viewed images counting novel+repeat presentations
end
valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

binsize = 12; %pixels per bin spatial bin in either dimension 1/2 dva
filter_width = 6; %std of 2D guassian filter ~ 3 dva, could use 2 dva (4) as well
H = define_spatial_filter(filter_width);

switch task
    case {'cvtnew','CVTNEW'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to dot position---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        info_type = 'spatial_cvtnew';%type of mutual information analysis to perform
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'meta');
        
        %get a generic calibration function for now
        if strcmpi(task_file(1:2),'TO')
            load([data_dir 'TO151203_3-preprocessed.mat'],'tform');
        else
            load([data_dir 'PW140729_3-preprocessed.mat'],'tform');
        end
        
        [eyechans] = find_desired_channels(cfg,'eye');

        
        %slight bug in the timing file that reporst the wrong condition
        %number for a few trials on a few sessions which was later fixed,
        %removes 1 trial/block or path
        if strcmpi(task_file(1:2),'PW)')
            file_date = str2double(task_file(3:6));
            if file_data >= 141007 && file_data <= 141024
                [meta,cfg] = CVTNEW_remove_bad_condition(meta,cfg);
            end
        end
        
        disp('Collecting Spike Locations')
        
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        dot_on_code = 25;
        dot_clrchng_code = 27;
        bar_code_response = 4; %monkey made a move
        reward_code = 3;
        imageX = 800;%528; %horizontal size of displayed dot positions
        imageY = 600;%528 %horizontal size of displayed dot positions
        early_response_code = 205;
        late_response_code = 202;
        break_fixation_code = 203;
        stimulus_off_code = 24;
        
        %preallocate space and parallel structure of cfg
        dotpos = cell(1,num_units);
        spike_times = cell(1,num_units);
        dotdirection = cell(1,num_units);
        correct_trials = cell(1,num_units);
        cross_eyepos = cell(2,num_units);
        for unit = 13%1:num_units
            spike_times{unit} = NaN(length(cfg.trl),3000);
            dotpos{unit} = NaN(2*length(cfg.trl),3000);
            dotdirection{unit} =  NaN(length(cfg.trl),3000);
            correct_trials{unit} = NaN(1,length(cfg.trl));
            cross_eyepos{1,unit} = zeros(imageY,imageX);
            cross_eyepos{2,unit} = zeros(imageY,imageX);
        end
        
        for t = 1:length(cfg.trl);
            if any(cfg.trl(t).allval == dot_on_code);% any(cfg.trl(t).allval == reward_code); %only take correct trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                pathon =  cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %meta data starts at 1st 666 which appears before dot actually turns on in event array
                dot_clrchng = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start;
                responded = cfg.trl(t).alltim(cfg.trl(t).allval == bar_code_response)-trial_start;
                correct = 1;
                if isempty(responded)%error trial
                    correct = 0;
                    responded = cfg.trl(t).alltim(cfg.trl(t).allval == early_response_code)-trial_start; %look for early response
                    if isempty(responded)
                        responded = cfg.trl(t).alltim(cfg.trl(t).allval == late_response_code)-trial_start; %look for early response
                        if isempty(responded)
                           responded = cfg.trl(t).alltim(cfg.trl(t).allval == break_fixation_code)-trial_start; %look for early response
                           if isempty(responded)
                                responded = cfg.trl(t).alltim(cfg.trl(t).allval == stimulus_off_code)-trial_start; %look for early response
                                responded = responded(1);
                           end
                        end
                    end
                end
                
                if length(meta(t).sample) < 50
                    continue
                end
                if length(meta(t).sample) < 250
                    xn = interp1(meta(t).sample,meta(t).x,meta(t).sample(1):meta(t).sample(end),'linear');
                    yn = interp1(meta(t).sample,meta(t).y,meta(t).sample(1):meta(t).sample(end),'linear');
                else
                    xn = interp1(meta(t).sample,meta(t).x,meta(t).sample(1):meta(t).sample(end),'cubic');
                    yn = interp1(meta(t).sample,meta(t).y,meta(t).sample(1):meta(t).sample(end),'cubic');
                end
                direction =  atan2(diff(yn),diff(xn));
                direction = [direction direction(end)]; %correlated in time and want to make equal size to spike train
                
                for unit = 13%1:num_units
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                        spikes = find(data(unit).values{t});
                        spikeind = spikes(spikes > pathon & spikes <= responded)-pathon;
                        spikeind(spikeind < 1) = []; %should only happen when spikes occur at the same time as tstart
                        
                        if correct
                            %generic calibration for now
                            xeye = data(eyechans(1)).values{t}(pathon:responded);
                            yeye = data(eyechans(2)).values{t}(pathon:responded);
                            [xeye,yeye] = tformfwd(tform,xeye,yeye); 
                            
                            %convert from dva to pixels
                            xeye = 24*xeye; 
                            xeye = xeye+imageX/2;
                            xeye = round(xeye);
                            yeye = 24*yeye;
                            yeye = yeye+imageY/2;
                            yeye = round(yeye);

                            for xy = 1:length(xeye)
                                cross_eyepos{1,unit}(yeye(xy),xeye(xy)) = cross_eyepos{1,unit}(yeye(xy),xeye(xy)) + 1; %all time
                            end
                            
                            %get eye pos when spike occurs
                            for s = 1:length(spikeind)
                                cross_eyepos{2,unit}(yeye(spikeind(s)),xeye(spikeind(s))) =...
                                    cross_eyepos{2,unit}(yeye(spikeind(s)),xeye(spikeind(s))) +1;%just when spikes occur
                            end
                        end
                        
                        if all(isnan(xn)) && correct == 0
                            %some error must of occured can't use this data
                            %so continue and not worth it either
                            continue
                        elseif all(isnan(xn)) && correct == 1
                            error('wtf')
                        end
                            
                                                    
                        correct_trials{unit}(t) = correct;
                        
                        temp = zeros(1,length(xn));
                        temp(spikeind) = 1;
                        spike_times{unit}(t,1:length(temp)) = temp;
                        
                        dotpos{unit}(2*t-1,1:length(xn)) = round(xn);
                        dotpos{unit}(2*t,  1:length(xn)) = round(yn);
                        dotdirection{unit}(t,1:length(xn)) = direction;
                    end
                end
            end
        end
        
        %remove excess NaNs associated with error trials
        dotpos = laundry(dotpos);
        correct_trials=laundry(correct_trials);
        spike_times = laundry(spike_times);
        dotdirection = laundry(dotdirection);
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Determining if neurons are spatially modulated')
        
        %---Determine if Neuron is more spatial than chance---%
        spatial_info.rate = NaN(1,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(1,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_rate_prctile = NaN(1,num_units);
        
        %spatial correlation first vs second half row 1 spearman row 2 Kendall
        spatial_info.spatialstability_halves = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_halves = cell(1,num_units);
        spatial_info.spatialstability_halves_prctile = NaN(1,num_units);
        
        %spatial correlation even and odd trials row 1 spearman row 2 Kendall
        spatial_info.spatialstability_even_odd = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_even_odd = cell(1,num_units);
        spatial_info.spatialstability_even_odd_prctile = NaN(1,num_units);
        
        trial_data{3} = [imageX imageY];
        for unit = 13%1:num_units
            if ~isempty(dotpos{unit})
                %remove data for incorrect trials since may not have been paying attention
                dtpsx = dotpos{unit}(1:2:end,:);
                dtpsx = dtpsx(correct_trials{unit} == 1,:);
                dtpsy = dotpos{unit}(2:2:end,:);
                dtpsy = dtpsy(correct_trials{unit} == 1,:);
                
                dtps = NaN(2*size(dtpsx,1),size(dtpsx,2));
                dtps(1:2:end,:) = dtpsx;
                dtps(2:2:end,:) = dtpsy;

                trial_data{1} = dtps;
                trial_data{2} = spike_times{unit}(correct_trials{unit} == 1,:);
                
                [observed_info_rate,shuffled_info_rate]...
                    = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
                
                %skaggs information score
                spatial_info.rate(unit) = observed_info_rate.skaggs;
                spatial_info.shuffled_info_rate{unit} = shuffled_info_rate.skaggs;
                spatial_info.shuffled_rate_prctile(unit) = 100*sum(...
                    spatial_info.rate(unit) > spatial_info.shuffled_info_rate{unit})/numshuffs;
                
                %spatial correlation first half vs second half
                spatial_info.spatialstability_halves(unit) = observed_info_rate.spatialstability_halves;
                spatial_info.shuffled_spatialstability_halves{unit} = shuffled_info_rate.spatialstability_halves;
                spatial_info.spatialstability_halves_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_halves(unit) > ...
                    spatial_info.shuffled_spatialstability_halves{unit})/numshuffs;
                
                %spatial correlation for even and odd trials
                spatial_info.spatialstability_even_odd(unit) = observed_info_rate.spatialstability_even_odd;
                spatial_info.shuffled_spatialstability_even_odd{unit} =...
                    shuffled_info_rate.spatialstability_even_odd;
                spatial_info.spatialstability_even_odd_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_even_odd(unit) > ...
                    spatial_info.shuffled_spatialstability_even_odd{unit})/numshuffs;
            end
        end
        
        path_numbers = cell(1,length(num_units));
        for unit = 1:num_units
            all_cnds = NaN(1,length(cfg.trl));
            for t = 1:length(cfg.trl);
                if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only put data in for valid trials
                    all_cnds(t) = cfg.trl(t).cnd-1000;
                end
            end
            all_cnds(isnan(all_cnds)) = [];
            unique_cnds = unique(all_cnds,'stable');
            paths = NaN(1,length(cfg.trl));
            for c = 1:length(unique_cnds)
                paths(all_cnds == unique_cnds(c)) = c;
            end
            path_numbers{unit} = paths; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----Determine if Neuron is Directionally Modulated---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shuffled_95_direction_STA = cell(1,num_units);
        observed_STA = cell(1,num_units);
        observed_direction = cell(1,num_units);
        for unit = 1:num_units
            [shuffled_95_direction_STA{unit},observed_STA{unit},observed_direction{unit}] = ...
                CVTNEW_Direction_Analysis(spike_times{unit},dotdirection{unit},numshuffs,correct_trials{unit});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task_type = 'cvtnew_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = unit_stats(1,:);
        unit_names.path_number = path_numbers;
        unit_names.correct_trials = correct_trials;
        unit_names.eyepos = cross_eyepos;
        spatial_info.shuffled_95_direction_STA = shuffled_95_direction_STA;
        spatial_info.observed_STA = observed_STA;
        spatial_info.observed_direction = observed_direction;
        spatial_analysis_plotsV2(figure_dir,task_file,dotpos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,NaN,Fs);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        save([data_dir task_file(1:10) '-spatial_analysis_results.mat'],...
            'spike_times','dotpos','spatial_info','binsize','filter_width',...
            'dotdirection','unit_names','shuffled_95_direction_STA','observed_STA',...
            'observed_direction')
        disp(['Spatial Data Analyis for ' task_file(1:10) ' saved']);
        
    case 'ListSQ'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Reformat data so that spikes are locked to eye movements---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        image_on_twin = 500;%how much time to ignore eye movements for to remove
        %strong visual response (though some may last longer) combined with strong central bias
        
        info_type = 'spatial';%type of mutual information analysis to perform
        
        disp('Collecting Spike Locations')
        
        %set/get some general important info
        [eyechans] = find_desired_channels(cfg,'eye');
        num_trials = length(cfg.trl);
        Fs = data(1).fsample; %should be 1000
        ITIstart_code = 15;
        img_on_code= 23;
        img_off_code = 24;
        imageX = 800; %horizontal size of the image in pixels
        imageY = 600; %horizontal size of the image in pixels
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805 %first sets done by PW before this data had 7 s images
            imgdur = 7000;
        else %rest of image presentations were 5 seconds
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        %---preallocate space and parallel structure of cfg---%
        %Each unit gets own cell since each unit may have different amount of trials/data
        eyepos = cell(1,num_units);
        spike_times = cell(1,num_units);
        which_images = cell(1,num_units);
        nvr = cell(1,num_units);
        for unit = 1:num_units
            spike_times{unit} = NaN(length(cfg.trl),imgdur*1.5);
            eyepos{unit} = NaN(2*length(cfg.trl),imgdur*1.5);
            which_images{unit} = NaN(1,length(cfg.trl));
            nvr{unit} = NaN(1,length(cfg.trl));
        end
        
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code); %should be first code in sequence
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %want to avoid visual response
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start;
                
                % if monkey isn't paying attention and looked away image presentation
                % is now longer than imgdur (because of cumulative looking time)
                % so data isn't probably worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                    imgoff = imgon+1.5*imgdur-1;
                end
                imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                
                xn = round(data(eyechans(1)).values{t}(imgon:imgoff)); %horizontal eye data
                yn = round(data(eyechans(2)).values{t}(imgon:imgoff)); %vertical eye data
                
                img_index = find(img_cnd == cfg.trl(t).cnd); %get image index
                
                if any(isnan(which_img(img_index)))%presentation error so skip trial, see get_image_numbers.m
                    continue
                end
                
                for unit = 1:num_units
                    if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                        eyepos{unit}(2*t-1,1:length(xn)) = xn; %horizontal eye data
                        eyepos{unit}(2*t,1:length(yn)) = yn; %vertical eye data
                        
                        which_images{unit}(t) = which_img(img_index);
                        nvr{unit}(t) = novel_vs_repeat(img_index);
                        
                        spikes = find(data(unit).values{t});
                        if ~isempty(spikes)
                            spikeind = spikes(spikes > imgon & spikes <= imgoff)-imgon;
                            temp = zeros(1,length(xn));
                            temp(spikeind) = 1;
                            spike_times{unit}(t,1:length(temp)) = temp;
                        else
                            temp = zeros(1,length(xn));
                            spike_times{unit}(t,1:length(temp)) = temp;
                        end
                    end
                end
            end
        end
        
        %remove excess NaNs associated with error trials
        eyepos = laundry(eyepos,1); %only remove rows aka trials without data
        spike_times = laundry(spike_times,1);  %only remove rows aka trials without data
        which_images = laundry(which_images);
        nvr = laundry(nvr);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Perform Spatial analysis and signifance testing---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Determining if neurons are spatially modulated')
        
        %--calculate peak firing rate---%
        peak_firing_rate = NaN(3,num_units);%row 1 novel, row 2 repeat, row 3 all images
        for unit = 1:num_units
            if isnan(valid_trials(1,unit)); %no valid trials for unit
                continue
            end
            for condition = 1:3
                if condition == 1 %novel images
                    trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 1);
                    trial_data{2} = spike_times{unit}(nvr{unit} == 1,:);
                elseif condition == 2 %repeat images
                    trial_data{1} = select_eyepos(eyepos{unit},nvr{unit} == 2);
                    trial_data{2} = spike_times{unit}(nvr{unit} == 2,:);
                elseif condition == 3 %all images
                    trial_data{1} = eyepos{unit};
                    trial_data{2} = spike_times{unit};
                end
                if nansum(nansum(trial_data{2})) == 0 %this occassiionally happens with really specifc/low firing rate neurons...
                    %and we don't want a NaN
                    peak_firing_rate(condition,unit) = 0;
                    continue
                end
                
                [ratemap,~] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,'all');
                peak_firing_rate(condition,unit) = prctile(ratemap(:),97.5);%don't want to grab outliers
            end
        end
        
        %---Determine if Neuron is more spatial than chance---%
        spatial_info.rate = NaN(1,num_units); %the observed information rate in bits/sec
        spatial_info.shuffled_info_rate = cell(1,num_units); %bootstrapped information rate in bits/sec expected by chance from spike train
        spatial_info.shuffled_rate_prctile = NaN(1,num_units);
        
        %spatial correlation first vs second half row 1 spearman row 2 Kendall
        spatial_info.spatialstability_halves = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_halves = cell(1,num_units);
        spatial_info.spatialstability_halves_prctile = NaN(1,num_units);
        
        %spatial correlation even and odd trials row 1 spearman row 2 Kendall
        spatial_info.spatialstability_even_odd = NaN(1,num_units);
        spatial_info.shuffled_spatialstability_even_odd = cell(1,num_units);
        spatial_info.spatialstability_even_odd_prctile = NaN(1,num_units);
        
        trial_data{3} = [imageX imageY];
        for unit = 1:num_units
            if ~isempty(spike_times{unit})
                trial_data{1} = eyepos{unit};
                trial_data{2} = spike_times{unit};
                [observed_info_rate,shuffled_info_rate]...
                    = estimated_mutual_information(trial_data,numshuffs,info_type,[binsize,filter_width],Fs);
                
                %skaggs information score
                spatial_info.rate(unit) = observed_info_rate.skaggs;
                spatial_info.shuffled_info_rate{unit} = shuffled_info_rate.skaggs;
                spatial_info.shuffled_rate_prctile(unit) = 100*sum(...
                    spatial_info.rate(unit) > spatial_info.shuffled_info_rate{unit})/numshuffs;
                
                %spatial correlation first half vs second half
                spatial_info.spatialstability_halves(unit) = observed_info_rate.spatialstability_halves;
                spatial_info.shuffled_spatialstability_halves{unit} = shuffled_info_rate.spatialstability_halves;
                spatial_info.spatialstability_halves_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_halves(unit) > ...
                    spatial_info.shuffled_spatialstability_halves{unit})/numshuffs;
                
                
                %spatial correlation for even and odd trials
                spatial_info.spatialstability_even_odd(unit) = observed_info_rate.spatialstability_even_odd;
                spatial_info.shuffled_spatialstability_even_odd{unit} =...
                    shuffled_info_rate.spatialstability_even_odd;
                spatial_info.spatialstability_even_odd_prctile(unit) = ...
                    100*sum(spatial_info.spatialstability_even_odd(unit) > ...
                    spatial_info.shuffled_spatialstability_even_odd{unit})/numshuffs;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot and Save Figures of Results---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        task_type = 'List_spatial';
        unit_names.multiunit = multiunit;
        unit_names.name = unit_stats(1,:);
        unit_names.avg_firing_rate = whole_session_mean_firing_rate;
        unit_names.putative_EI = excitatory_inhibitory;
        spatial_analysis_plotsV2(figure_dir,task_file,eyepos,spike_times,spatial_info,task_type,...
            unit_names,[binsize,filter_width],imageX,imageY,nvr,Fs,peak_firing_rate);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Finally save all the data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],...
            'spike_times','eyepos','spatial_info','nvr','which_images',...
            'binsize','filter_width','peak_firing_rate','image_on_twin','Fs',...
            'whole_session_mean_firing_rate','excitatory_inhibitory','unit_names')
        disp(['ListSQ Spatial Data Analyis for ' task_file(1:end-11) ' saved']);
end
end