%Code Written to Make Time Cell Plots Limited by Fixation Duration Rather
%than blindly plotting from fixation start to twin (usually 400 ms).
%Code Is essentially the same as Population_Place_Cell_Fixation_Anlaysis
%otherwise.
%Written by Seth Koenig 5/15/19

clar %clear,clc

task = 'ListSQ';
min_blks = 2; %only anaedlyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial

%---Misc. Parameters (Mostly Place Cells)---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_place_cell_unit_names = {}; %place cell unit names
all_place_cell_monkeys = []; %1s and 2s

%---Spatial Correlation Values for all cells---%
all_place_cell_spatial_corrs = []; %place cell spatial correlations
all_place_cell_spatial_skagg_percents= []; %place cell skagg information score percentiles
all_non_place_cell_spatial_corrs = [];%non-place cell spatial correlations
all_non_place_cell_skagg_percents = [];%non-place cell skagg information score percentiles
non_place_skagg_names = [];

%---Place Cell Firing Rate Curves for List---%
all_list_peak_times = [];%all peak times
all_out2in_rates = [];%fixations out-> in for -twin:2*twin
all_out2in_rates_limited = [];%fixations out-> in limited by fixation duration
all_out2in_rates_limited2 = [];%fixations out-> in limited by fixation duration
all_out2in_rates_limited3 = [];%limited2 plus removoval of short fix before out
all_out2in_minus_out2out_rates = [];
all_out2in_minus_out2out_rates2 = [];

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
        
        %load spatial analysis data
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
        
        for unit = 1:num_units
            if ~isempty(spatial_info.shuffled_info_rate{unit})&& length(spatial_info.shuffled_info_rate{unit}) < 1000
                error('Should have 1000 shuffles') %make sure file has right number of shuffles
            elseif isempty(spatial_info.shuffled_info_rate{unit})
                continue
            end
            
            %             if isnan(stats_across_tasks(1,unit))%no peak detected
            %                 continue
            %             end
            %
            if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
                    && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
                
                
                %---Misc Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_place_cell_unit_names = [all_place_cell_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_place_cell_monkeys = [all_place_cell_monkeys monk]; %1s and 2s for monkey
                
                %---Spatial Correlation Values for Place cells---%
                all_place_cell_spatial_corrs = [all_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_place_cell_spatial_skagg_percents = [all_place_cell_spatial_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                
                
                place_field_matrix = all_place_field_matrix{unit};
                %note ignores incomplete coverage

                %---Out2In calculated from original source (as a double-check)---%
                %firing rate out-> in
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 1,:); %get spike trains
                out2in_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                
                %firing rate out-> out
                firing_rate = list_fixation_locked_firing{unit}(in_out{unit} == 4,:); %get spike trains
                out2out_curve = nandens(firing_rate,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                out2in_minus_out2out_curve = out2in_curve-out2out_curve;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%---Calculate Firing Rate Locked to Fixations---%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                fix_in_out = NaN(1,1000); %firing rate for fixations inside or outside  of place field
                %1) first fixation in: out-> in
                %2) fixation in field but not first: in -> in
                %3) first fixation out of field: in -> out
                %4) fixation out of field but not first: out-> out
                
                next_fix_in_out = NaN(1,1000);
                fix_locked_firing = NaN(1000,twin1+twin2); %spike trains locked to fixations
                all_fix_dur = NaN(1,1000);
                all_previous_fix_dur = NaN(1,1000);%include prior sacdur
                
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
                                
                                %determine if next fixation was in or out of place field
                                if f~= size(fixations,2)
                                    fixx = fixations(1,f+1);
                                    fixx(fixx < 1) = 1;
                                    fixx(fixx > imageX) = imageX;
                                    fixy = imageY-fixations(2,f+1);
                                    fixy(fixy < 1) = 1;
                                    fixy(fixy > imageY) = imageY;
                                    
                                    if place_field_matrix(fixy,fixx) == 1 %then prior fixation was already inside
                                        next_fix_in_out(fix_ind) = 1;
                                    else
                                        next_fix_in_out(fix_ind) = 0;
                                    end
                                    last_fixx = fixx;
                                    last_fixy = fixy;
                                end
                                
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
                                
                                if fix_in_out(fix_ind) == 1 %out2in & next fix is not in field
                                    all_fix_dur(fix_ind) = fixdur;
                                    all_previous_fix_dur(fix_ind) = fixationtimes(1,f)-fixationtimes(1,f-1);
                                    %get firing rate locked to fixation
                                    fixt = fixationtimes(1,f);%start of fixation
                                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                                    temp = zeros(1,twin1+twin2);
                                    temp(fix_spikes) = 1;
                                    fix_locked_firing(fix_ind,:) = temp;
                                    
                                    fix_ind = fix_ind+1;
                                end
                            end
                        end
                    end
                end
                %%
                [out2in_curve_recomp,out_2in_smoothed_trials] = nandens(fix_locked_firing,smval,'gauss',Fs,'nanflt'); %calculate smoothed firing rate
                
                out2in_curve_limited = out_2in_smoothed_trials;
                out2in_curve_limited2 = fix_locked_firing;
                out2in_curve_limited3 = fix_locked_firing;
                for fix = 1:fix_ind
                    if all_fix_dur(fix) < twin2 && next_fix_in_out(fix) == 1
                        out2in_curve_limited(fix,twin1+all_fix_dur(fix)+1:end) = NaN;%remove time by start of next fixation
                        out2in_curve_limited2(fix,twin1+all_fix_dur(fix)+1:end) = NaN;%remove time by start of next fixation
                        if all_fix_dur(fix)+44 < twin2
                            out2in_curve_limited3(fix,twin1+all_fix_dur(fix)+44:end) = NaN;%remove time by start of next fixation
                        end
                    end
                end
                n_trials = sum(~isnan(out2in_curve_limited));
                out2in_curve_limited = nandens(nanmean(out2in_curve_limited),smval/3,'gauss',Fs,'nanflt')/1000;
                n_trials_too_little = find(n_trials < 50);
                %out2in_curve_limited(n_trials_too_little) = NaN;
                
                out2in_curve_limited2 = nandens3(out2in_curve_limited2,smval,Fs); %calculate smoothed firing rate
                out2in_curve_limited3 = nandens3(out2in_curve_limited3,smval,Fs); %calculate smoothed firing rate
                %%
%                 figure
%                 plot(-twin1:twin2-1,out2in_curve)
%                 hold on
%                 plot(-twin1:twin2-1,out2out_curve)
%                 plot(-twin1:twin2-1,out2in_curve_recomp-0.05);
%                 plot(-twin1:twin2-1,out2in_curve_limited);
%                 plot(-twin1:twin2-1,out2in_curve_limited2);
%                 plot(-twin1:twin2-1,out2in_curve_limited3);
%                 hold off
%                 xlabel('Time from Fixation Start (ms)')
%                 ylabel('Firing Rate')
%                 legend({'Out2In','Out2Out','Recomputed Out2In','Smooth then Avg','Avg Then Smooth','Avg Then Smooth 3'},...
%                     'location','NorthEastOutside')
                %%
                
                %---Normalize---%
                out2in_curve_recomp = out2in_curve_recomp-nanmean(out2in_curve_recomp(1:twin1)); %remove base line
                out2in_curve_recomp = out2in_curve_recomp/max(out2in_curve_recomp);
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    out2in_curve_recomp = out2in_curve_recomp/out2in_curve_recomp(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    out2in_curve_recomp = out2in_curve_recomp/max(out2in_curve_recomp);
                end
                
                out2in_curve_limited = out2in_curve_limited-nanmean(out2in_curve_limited(1:twin1)); %remove base line
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    out2in_curve_limited = out2in_curve_limited/out2in_curve_limited(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    out2in_curve_limited = out2in_curve_limited/max(out2in_curve_limited);
                end
                
                out2in_curve_limited2 = out2in_curve_limited2-nanmean(out2in_curve_limited2(1:twin1)); %remove base line
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    out2in_curve_limited2 = out2in_curve_limited2/out2in_curve_limited2(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    out2in_curve_limited2 = out2in_curve_limited2/max(out2in_curve_limited2);
                end
                
                out2in_curve_limited3 = out2in_curve_limited3-nanmean(out2in_curve_limited3(1:twin1)); %remove base line
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                    out2in_curve_limited3 = out2in_curve_limited3/out2in_curve_limited3(stats_across_tasks(1,unit));%divide by peak firing rate
                else %no peak so normalize by max
                    out2in_curve_limited3 = out2in_curve_limited3/max(out2in_curve_limited3);
                end
                
                out2in_minus_out2out_curve2 = out2in_minus_out2out_curve-nanmean(out2in_minus_out2out_curve(1:twin1));
                out2in_minus_out2out_curve3 = out2in_minus_out2out_curve2;
                if ~isnan(stats_across_tasks(1,unit))%peak detected
                     out2in_minus_out2out_curve2 = out2in_minus_out2out_curve2/out2in_minus_out2out_curve2(stats_across_tasks(1,unit));%divide by peak firing rate
                else
                    out2in_minus_out2out_curve2 = out2in_minus_out2out_curve2/max(out2in_minus_out2out_curve2);
                end
                
                out2in_minus_out2out_curve3 = out2in_minus_out2out_curve3/max(out2in_minus_out2out_curve3);
                
                all_list_peak_times = [all_list_peak_times stats_across_tasks(1,unit)];%time of peak firing rate of place cells for out-> in fixations
                all_out2in_rates = [all_out2in_rates; out2in_curve_recomp];%fixations out-> in for -twin:2*twin
                all_out2in_rates_limited = [all_out2in_rates_limited; out2in_curve_limited];%fixations out-> in limited by fixation duration
                all_out2in_rates_limited2 = [all_out2in_rates_limited2; out2in_curve_limited2];%fixations out-> in limited by fixation duration
                all_out2in_rates_limited3 = [all_out2in_rates_limited3; out2in_curve_limited3];%fixations out-> in limited by fixation duration
                all_out2in_minus_out2out_rates = [all_out2in_minus_out2out_rates; out2in_minus_out2out_curve2];
                all_out2in_minus_out2out_rates2 = [all_out2in_minus_out2out_rates2; out2in_minus_out2out_curve3];

            elseif ~isnan(spatial_info.shuffled_rate_prctile(unit))
                
                %---Misc Parameters+---%
                monkey_all_unit_count(2,monk) = monkey_all_unit_count(2,monk)+1; %unit count
                all_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                
                %---Spatial Correlation Values for Non-Place cells---%
                all_non_place_cell_spatial_corrs = [all_non_place_cell_spatial_corrs spatial_info.spatialstability_halves(unit)];
                all_non_place_cell_skagg_percents = [all_non_place_cell_skagg_percents spatial_info.shuffled_rate_prctile(unit)];
                non_place_skagg_names = [non_place_skagg_names  {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                
            end
        end
    end
end
%%
%remove neurons without a definitive peak
all_out2in_rates(isnan(all_list_peak_times),:) = [];
all_out2in_rates_limited(isnan(all_list_peak_times),:) = [];
all_out2in_rates_limited2(isnan(all_list_peak_times),:) = [];
all_out2in_rates_limited3(isnan(all_list_peak_times),:) = [];
all_list_peak_times(isnan(all_list_peak_times)) = [];
%%
[~,place_order] = sort(all_list_peak_times);

% [mx,mxi] = max(all_out2in_rates');
% [~,place_order] = sort(mxi);

figure
subplot(2,2,1)
vals = all_out2in_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_rates,1)],all_out2in_rates(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Original: out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square

subplot(2,2,2)
vals = all_out2in_rates_limited(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_rates_limited,1)],all_out2in_rates_limited(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_rates_limited,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Limited Smooth then Avg: out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square

subplot(2,2,3)
vals = all_out2in_rates_limited2(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_rates_limited2,1)],all_out2in_rates_limited2(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_rates_limited2,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Limited Avg then Smooth2: out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square

subplot(2,2,4)
vals = all_out2in_rates_limited3(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_rates_limited3,1)],all_out2in_rates_limited3(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_rates_limited3,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Limited Avg then Smooth3: out -> in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
axis square
%%
figure
subplot(2,2,1)
vals = all_out2in_minus_out2out_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_minus_out2out_rates,1)],all_out2in_minus_out2out_rates(place_order,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_minus_out2out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2Out - Out2in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

[~,mxi] = max(all_out2in_minus_out2out_rates');
[~,place_order2] = sort(mxi);
subplot(2,2,3)
vals = all_out2in_minus_out2out_rates(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_minus_out2out_rates,1)],all_out2in_minus_out2out_rates(place_order2,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_minus_out2out_rates,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2Out - Out2in')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar

[~,mxi] = max(all_out2in_minus_out2out_rates2');
[~,place_order3] = sort(mxi);
subplot(2,2,2)
vals = all_out2in_minus_out2out_rates2(:,1:twin1); %"baseline" out of field firing rate
imagesc([-twin1:twin2-1],[1:size(all_out2in_minus_out2out_rates2,1)],all_out2in_minus_out2out_rates2(place_order3,:))
colormap('jet')
hold on
plot([0 0],[1 size(all_out2in_minus_out2out_rates2,1)],'w--');
hold off
xlabel('Time from Fixation Start')
ylabel('Neuron #')
title('Out2Out - Out2in Self-Norm')
caxis([-std(vals(:)) 1]) %set minumun to standard deviation of baseline since some neurons are greatly inhibited
xlim([-twin1 twin2])
colorbar
%%

figure
hold on
plot([-twin1:twin2-1],mean(all_out2in_rates))
plot([-twin1:twin2-1],mean(all_out2in_rates_limited))
plot([-twin1:twin2-1],mean(all_out2in_rates_limited2))
plot([-twin1:twin2-1],mean(all_out2in_rates_limited3))
plot([-twin1:twin2-1],mean(all_out2in_minus_out2out_rates))
plot([-twin1:twin2-1],mean(all_out2in_minus_out2out_rates2))
plot([-twin1 twin2-1],[0 0],'k')
plot([0 0],[-0.2 1],'k--')
plot([-44 -44],[-0.2 1],'k--')
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Normalized Firing Rate')
legend('Original','Limited SMooth then Avg','Limited Avg then smooth','Limited Avg Then Smooth3','Out2in - Out2out','Out2in - Out2out Self-norm')
%%
