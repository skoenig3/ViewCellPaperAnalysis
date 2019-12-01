% Code below creates population summary for Significnat Place Cells
% Written by Seth Konig
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

%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Desktop\Significant Units\Spatial Analysis\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

imageX = 800;
imageY = 600;
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial
numshuffs = 10000; %number of shuffles for resmampling

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

Fs = 1000; %Hz sampling frequency
fixwin = 5;%size of fixation window on each crosshair
smval = 30;%2*std of gaussian kernel so 15 ms standard deviation
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer


%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey

all_peak_list = [];
%---Place Cell Firing Rate Curves for List---%
all_in_rates = [];%fixations out-> in
all_in_rates100 = [];%fixations out-> in
all_in_rates200 = [];%fixations out-> in
all_in_rates300 = [];%fixations out-> in

all_in_novel = [];
all_in_repeat = [];

sig_nov_rep_out2in = [];
sig_nov_rep_stimulus = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk =2:-1:1
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
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
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
        
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']) %visual response analysis
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results']) %memory visual response analysis
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        %load Place Cell Fixation Analysis data
        load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])
        if numshuffs < 1000
            error('Should have 1000 shuffles')%make sure file has right number of shuffles
        end
        if smval ~= 30
            error('Smoothing Value (2xStd) does not match expectations!')
        end
        
        if any(spatial_info.shuffled_rate_prctile> 95  & spatial_info.spatialstability_halves_prctile > 95)
        else
            continue
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
                
                place_field_matrix = all_place_field_matrix{unit};
                %note ignores incomplete coverage
                
                %---for image on short 1 second window---%
                if (epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)
                    sig_nov_rep_stimulus = [sig_nov_rep_stimulus 2];
                elseif (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    sig_nov_rep_stimulus = [sig_nov_rep_stimulus 1];
                else
                    sig_nov_rep_stimulus = [sig_nov_rep_stimulus 0];
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Calculate Firing Rate Locked to Fixations---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
                %1) first fixation in: out-> in
                %2) fixation in field but not first: in -> in
                %3) first fixation out of field: in -> out
                %4) fixation out of field but not first: out-> out
                fix_locked_firing = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                fix_locked_firing100 = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                fix_locked_firing200 = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                fix_locked_firing300 = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                
                fix_locked_firing_nov = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                fix_locked_firing_rep = NaN(1000,(twin1+twin2)); %spike trains locked to fixations
                
                fix_nov_ind = 1;
                fix_rep_ind = 1;
                
                fix_ind = 1; %fixation # so can track in variables above
                fix_ind100 = 1; %fixation # so can track in variables above
                fix_ind200= 1; %fixation # so can track in variables above
                fix_ind300 = 1; %fixation # so can track in variables above
                
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
                            novrep= novel_vs_repeat(img_index); %whether image was novel or repeated
                            
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
                                        continue
                                    else %out->in
                                        fix_in_out(fix_ind) = 1;
                                    end
                                elseif place_field_matrix(fixy,fixx) == 0 %then inside, NaNs are for locations not analyzed
                                    if prior_fix_in_out == 1%prior fixation was inside so in->out
                                        fix_in_out(fix_ind) = 3;
                                        continue
                                    else %out->out
                                        fix_in_out(fix_ind) = 4;
                                        continue
                                    end
                                else %not a valid fixation location too sparse of coverage to analyze
                                    continue
                                end
                                
                                %get firing rate locked to fixation
                                fixt = fixationtimes(1,f);%start of fixation
                                fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                                temp = zeros(1,twin1+twin2);
                                temp(fix_spikes) = 1;
                                fix_locked_firing(fix_ind,:) = temp;
                                fix_ind = fix_ind+1;
                                
                                if novrep == 1
                                    fix_locked_firing_nov(fix_nov_ind,:) = temp;
                                    fix_nov_ind = fix_nov_ind+1;
                                else
                                    fix_locked_firing_rep(fix_rep_ind,:) = temp;
                                    fix_rep_ind = fix_rep_ind+1;
                                end
                                
                                if fixdur < 168 %automatically greater than 100
                                    fix_locked_firing100(fix_ind100,:) = temp;
                                    fix_ind100 = fix_ind100+1;
                                elseif fixdur < 222 %& > 200
                                    fix_locked_firing200(fix_ind200,:) = temp;
                                    fix_ind200 = fix_ind200+1;
                                else
                                    fix_locked_firing300(fix_ind300,:) = temp;
                                    fix_ind300 = fix_ind300+1;
                                end
                            end
                        end
                    end
                end
                
                %remove excess NaNs
                fix_locked_firing = laundry(fix_locked_firing);
                fix_locked_firing100 = laundry(fix_locked_firing100);
                fix_locked_firing200 = laundry(fix_locked_firing200);
                fix_locked_firing300 = laundry(fix_locked_firing300);
                
                fix_locked_firing_nov =laundry(fix_locked_firing_nov);
                fr = nandens(fix_locked_firing_nov,smval,'gauss',Fs,'nanflt');
                fr = fr-nanmean(fr(1:twin1));
                fr = fr/nanmax(fr);
                all_in_novel = [all_in_novel; fr];
                
                fix_locked_firing_rep =laundry(fix_locked_firing_rep);
                fr = nandens(fix_locked_firing_rep,smval,'gauss',Fs,'nanflt');
                fr = fr-nanmean(fr(1:twin1));
                fr = fr/nanmax(fr);
                all_in_repeat = [all_in_repeat; fr];
                
                if size(fix_locked_firing_rep,1)+ size(fix_locked_firing_nov,1) ~= size(fix_locked_firing,1)
                    error('math does not compute')
                end
                
                %for fixations ?->in vs ?->out in case not enough samples for the one above
                novrep = [ones(1,size(fix_locked_firing_nov,1)) 2*ones(1,size(fix_locked_firing_rep,1))];
                [~,all_firing_rate_curve] = nandens([fix_locked_firing_nov; fix_locked_firing_rep],smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
                all_curves2 = NaN(numshuffs,twin1+twin2);
                
                parfor shuff = 1:numshuffs;
                    %for fixations out->in vs out->out
                    ind = randperm(length(novrep));
                    shuff_novrep = novrep(ind);%randomly assign fixations as in or out
                    
                    shuff_nov_curve = nanmean(all_firing_rate_curve(shuff_novrep == 1,:)); %out->in
                    shuff_rep_curve = nanmean(all_firing_rate_curve(shuff_novrep == 2,:));%out->out
                    all_curves(shuff,:) = shuff_nov_curve-shuff_rep_curve; %calculate shuffled difference
                end
                
                %for fixations out->in vs out->out
                nov_curve = nanmean(all_firing_rate_curve(novrep == 1,:));  %firing rate curve for fixations in the field
                rep_curve = nanmean(all_firing_rate_curve(novrep == 2,:));  %firing rate curve for fixations in the field
                observed_diff = nov_curve-rep_curve; %observed difference in firing rate
                [~,sig_times] = cluster_level_statistic(observed_diff,all_curves,1,smval); %multiple comparision corrected significant indeces
                
                if sum(sig_times) > 0
                    sig_nov_rep_out2in = [sig_nov_rep_out2in 1];
                    %%
                    figure
                    t12 = -200:1000-1;
                    t14 = -twin1:twin4-1;
                    
                    subplot(2,2,1)
                    hold on
                    dofill(t12,time_lock_firing{unit,3},'black',1,smval); %all trials
                    dofill(t12,time_lock_firing{unit,3}(nvr{unit} == 1,:),'blue',1,smval); %novel trials
                    dofill(t12,time_lock_firing{unit,3}(nvr{unit} == 2,:),'red',1,smval); %repeat trials
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Image On (ms)')
                    xlim([-twin1 twin2])
                    if epoch_data.rate_prctile(unit,3) > 95 && epoch_data.temporalstability_prctile(unit,3) > 95
                        title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
                            '\rho = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
                            num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
                    end
                    ylims(1,:) = ylim;
                    box off
                    
                    subplot(2,2,3)
                    hold on
                    dofill(t14,time_lock_firing{unit,5},'black',1,smval2); %all trials
                    dofill(t14,time_lock_firing{unit,5}(nvr{unit} == 1,:),'blue',1,smval2); %novel trials
                    dofill(t14,time_lock_firing{unit,5}(nvr{unit} == 2,:),'red',1,smval2); %repeat trials
                    hold off
                    ylabel('Firing Rate (Hz)')
                    xlabel('Time from Image On (ms)')
                    xlim([-twin1 twin4])
                    if epoch_data.rate_prctile(unit,5) > 95 && epoch_data.temporalstability_prctile(unit,5) > 95
                        title(['bit ' num2str(epoch_data.rate_prctile(unit,5),3) '% ' ...
                            '\rho = ' num2str(epoch_data.temporalstability(unit,5),2) ' ' ...
                            num2str(epoch_data.temporalstability_prctile(unit,5),3)  '%' ])
                    end
                    ylims(2,:) = ylim;
                    box off
                    subplot(2,2,2)
                    [trial,time] = find(fix_locked_firing_nov == 1);
                    [trial2,time2] = find(fix_locked_firing_rep == 1);
                    if ~isempty(trial)
                        trial2 = trial2 + max(trial);
                    end
                    plot(time-twin1,trial,'b.')
                    hold on
                    plot(time2-twin1,trial2,'r.')
                    plot([0 0],[0 max(trial2)],'k--')
                    hold off
                    xlim([-twin1 twin2])
                    ylim([0 max(trial2)])
                    box off
                    xlabel('Time From Fixation Start (ms)')
                    ylabel('Firing Rate (Hz)')
                    
                    subplot(2,2,4)
                    hold on
                    dofill(-twin1:twin2-1,fix_locked_firing_nov,'blue',1,smval);
                    dofill(-twin1:twin2-1,fix_locked_firing_rep,'red',1,smval);
                    yl = ylim;
                    if yl(1) < 0
                        yl(1) = 0;
                        ylim(yl);
                    end
                    plot([0 0],[yl(1) yl(2)],'k--')
                    gaps = findgaps(find(sig_times));
                    if ~isempty(gaps)
                        for g = 1:size(gaps,1)
                            gp = gaps(g,:);
                            gp(gp == 0) = [];
                            h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                                [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                            uistack(h,'down')
                            set(h,'facealpha',.25,'EdgeColor','None')
                        end
                    end
                    hold off
                    xlim([-twin1 twin2]);
                    xlabel('Time from Fixation Start (ms)')
                    ylabel('Firing Rate (Hz)')
                    legend('Novel','Repeat','Location','NorthWest')
                    subtitle([task_file ' ' unit_stats{1,unit}])
                    close all
                else
                    sig_nov_rep_out2in = [sig_nov_rep_out2in 0];
                end
                %%
                %                 figure
                %                 t12 = -200:1000-1;
                %                 t14 = -twin1:twin4-1;
                %
                %                 subplot(2,2,1)
                %                 hold on
                %                 dofill(t12,time_lock_firing{unit,3},'black',1,smval); %all trials
                %                 dofill(t12,time_lock_firing{unit,3}(nvr{unit} == 1,:),'blue',1,smval); %novel trials
                %                 dofill(t12,time_lock_firing{unit,3}(nvr{unit} == 2,:),'red',1,smval); %repeat trials
                %                 hold off
                %                 ylabel('Firing Rate (Hz)')
                %                 xlabel('Time from Image On (ms)')
                %                 xlim([-twin1 twin2])
                %                 if epoch_data.rate_prctile(unit,3) > 95 && epoch_data.temporalstability_prctile(unit,3) > 95
                %                     title(['bit ' num2str(epoch_data.rate_prctile(unit,3),3) '% ' ...
                %                         '\rho = ' num2str(epoch_data.temporalstability(unit,3),2) ' ' ...
                %                         num2str(epoch_data.temporalstability_prctile(unit,3),3)  '%' ])
                %                 end
                %                 ylims(1,:) = ylim;
                %                 box off
                %
                %                 subplot(2,2,3)
                %                 hold on
                %                 dofill(t14,time_lock_firing{unit,5},'black',1,smval2); %all trials
                %                 dofill(t14,time_lock_firing{unit,5}(nvr{unit} == 1,:),'blue',1,smval2); %novel trials
                %                 dofill(t14,time_lock_firing{unit,5}(nvr{unit} == 2,:),'red',1,smval2); %repeat trials
                %                 hold off
                %                 ylabel('Firing Rate (Hz)')
                %                 xlabel('Time from Image On (ms)')
                %                 xlim([-twin1 twin4])
                %                 if epoch_data.rate_prctile(unit,5) > 95 && epoch_data.temporalstability_prctile(unit,5) > 95
                %                     title(['bit ' num2str(epoch_data.rate_prctile(unit,5),3) '% ' ...
                %                         '\rho = ' num2str(epoch_data.temporalstability(unit,5),2) ' ' ...
                %                         num2str(epoch_data.temporalstability_prctile(unit,5),3)  '%' ])
                %                 end
                %                 ylims(2,:) = ylim;
                %                 box off
                %                 subplot(2,2,2)
                %                 [trial,time] = find(fix_locked_firing_nov == 1);
                %                 [trial2,time2] = find(fix_locked_firing_rep == 1);
                %                 if ~isempty(trial)
                %                     trial2 = trial2 + max(trial);
                %                 end
                %                 plot(time-twin1,trial,'b.')
                %                 hold on
                %                 plot(time2-twin1,trial2,'r.')
                %                 plot([0 0],[0 max(trial2)],'k--')
                %                 hold off
                %                 xlim([-twin1 twin2])
                %                 ylim([0 max(trial2)])
                %                 box off
                %                 xlabel('Time From Fixation Start (ms)')
                %                 ylabel('Firing Rate (Hz)')
                %
                %                 subplot(2,2,4)
                %                 hold on
                %                 dofill(-twin1:twin2-1,fix_locked_firing_nov,'blue',1,smval);
                %                 dofill(-twin1:twin2-1,fix_locked_firing_rep,'red',1,smval);
                %                 yl = ylim;
                %                 if yl(1) < 0
                %                     yl(1) = 0;
                %                     ylim(yl);
                %                 end
                %                 plot([0 0],[yl(1) yl(2)],'k--')
                %                 gaps = findgaps(find(sig_times));
                %                 if ~isempty(gaps)
                %                     for g = 1:size(gaps,1)
                %                         gp = gaps(g,:);
                %                         gp(gp == 0) = [];
                %                         h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                %                             [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
                %                         uistack(h,'down')
                %                         set(h,'facealpha',.25,'EdgeColor','None')
                %                     end
                %                 end
                %                 hold off
                %                 xlim([-twin1 twin2]);
                %                 xlabel('Time from Fixation Start (ms)')
                %                 ylabel('Firing Rate (Hz)')
                %                 legend('Novel','Repeat','Location','NorthWest')
                %                 subtitle([task_file ' ' unit_stats{1,unit}])
                %
                %                 save_and_close_fig('C:\Users\seth.koenig\Desktop\New folder\',[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
                
                %%
                
                all_peak_list = [all_peak_list stats_across_tasks(1,unit)];
                
                fr = nandens(fix_locked_firing,smval,'gauss',Fs,'nanflt');
                baseline = nanmean(fr(1:twin1));
                fr = fr-baseline;
                max_rate = fr(stats_across_tasks(1,unit));%divide by peak firing rate;
                fr = fr/max_rate;
                all_in_rates = [all_in_rates; fr]; %calculate smoothed firing rate];%fixations out-> in
                
                fr = nandens(fix_locked_firing100,smval,'gauss',Fs,'nanflt');
                %                 fr = fr-baseline;
                %                 fr = fr/max_rate;
                fr = fr-nanmean(fr(1:twin1));
                fr = fr/nanmax(fr);
                all_in_rates100 = [all_in_rates100; fr]; %calculate smoothed firing rate];%fixations out-> in
                
                fr = nandens(fix_locked_firing200,smval,'gauss',Fs,'nanflt');
                %                 fr = fr-baseline;
                %                 fr = fr/max_rate;
                fr = fr-nanmean(fr(1:twin1));
                fr = fr/nanmax(fr);
                all_in_rates200 = [all_in_rates200; fr]; %calculate smoothed firing rate];%fixations out-> in
                
                fr = nandens(fix_locked_firing300,smval,'gauss',Fs,'nanflt');
                %                 fr = fr-baseline;
                %                 fr = fr/max_rate;
                fr = fr-nanmean(fr(1:twin1));
                fr = fr/nanmax(fr);
                all_in_rates300 = [all_in_rates300; fr]; %calculate smoothed firing rate];%fixations out-> in
            else
                continue
            end
        end
    end
end
%%
figure
plot([-twin1:twin2-1],nanmean(all_in_rates))
hold on
plot([-twin1:twin2-1],nanmean(all_in_rates100))
plot([-twin1:twin2-1],nanmean(all_in_rates200))
plot([-twin1:twin2-1],nanmean(all_in_rates300))
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
plot([-44 -44],[yl(1) yl(2)],'k--')
plot([-twin1 twin2],[0 0],'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title('Average View Cell Fixation Aligned Firing Rate For Variable Fixation Durations')
legend('All','100-168 ms','168-222 ms','222+ ms')

%%
figure
subplot(2,2,1)
%[~,i] = max(all_in_rates,[],2); %sort by time of maximum firing rate
[~,all_order] = sort(all_peak_list);
vals = all_in_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_rates,1)],all_in_rates(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('All Out2In')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

subplot(2,2,2)
vals = all_in_rates100(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_rates100,1)],all_in_rates100(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates100,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2In 100-168 ms')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

subplot(2,2,3)
vals = all_in_rates200(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_rates200,1)],all_in_rates200(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates200,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2In 168-222 ms')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

subplot(2,2,4)
vals = all_in_rates300(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_rates300,1)],all_in_rates300(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_rates300,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2In 222+ ms')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar
%%

figure
subplot(2,2,1)
%[~,i] = max(all_in_rates,[],2); %sort by time of maximum firing rate;
vals = all_in_novel(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_novel,1)],all_in_novel(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_novel,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Novel Out2IN')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

subplot(2,2,3)
%[~,i] = max(all_in_rates,[],2); %sort by time of maximum firing rate;
vals = all_in_repeat(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_in_repeat,1)],all_in_repeat(all_order,:)) %sorted same order as above
colormap('jet')
hold on
plot([0 0],[1 size(all_in_repeat,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Repeat Out2IN')
caxis([-nanstd(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
colorbar

subplot(2,2,2)
plot([-twin1:twin2-1],nanmean(all_in_rates),'g')
hold on
plot([-twin1:twin2-1],nanmean(all_in_novel),'b')
plot([-twin1:twin2-1],nanmean(all_in_repeat),'r')
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--');
plot([-44 -44],[yl(1) yl(2)],'k--')
plot([-twin1 twin2],[0 0],'k')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
title('Average View Cell Population for Novel vs. Repeat')
legend('All','Novel','Repeat')

%%
prop = [sum((sig_nov_rep_stimulus == 1 | sig_nov_rep_stimulus == 2) & sig_nov_rep_out2in == 1)/length(sig_nov_rep_out2in) ...
    sum((sig_nov_rep_stimulus == 1 | sig_nov_rep_stimulus == 2) & sig_nov_rep_out2in == 0)/length(sig_nov_rep_out2in) ...
    sum(sig_nov_rep_stimulus == 0 & sig_nov_rep_out2in == 1)/length(sig_nov_rep_out2in)];

figure
subplot(2,2,1)
bar(1:3,100*prop)
hold on
plot([0 4],[5 5],'k--')
hold off
set(gca,'Xtick',[1 2 3])
set(gca,'XtickLabel',{'Stim & Field','Stim Only','Field Only'})
box off
ylabel('% of View Cells')
title('Novelty Effects to Stimulus Onset and View Field')

subplot(2,2,2)
bar(1:2,100*[sum(sig_nov_rep_stimulus == 1) sum(sig_nov_rep_stimulus == 2)]/length(sig_nov_rep_stimulus))
hold on
plot([0 3],[5 5],'k--')
hold off
set(gca,'XtickLabel',{'Short Window & Field','Long Window & Field'})
box off
ylabel('% of View Cells')
title('All View Cells: Short (1 sec) vs Long (5 sec) Window')

subplot(2,2,3)
bar(1:2,100*[sum(sig_nov_rep_stimulus == 1 & sig_nov_rep_out2in == 1) sum(sig_nov_rep_stimulus == 2 & sig_nov_rep_out2in == 1)]/length(sig_nov_rep_stimulus))
hold on
plot([0 3],[5 5],'k--')
hold off
set(gca,'XtickLabel',{'Short Window','Long Window'})
box off
ylabel('% of View Cells')
title('All View Cells: Short/Long & Field')