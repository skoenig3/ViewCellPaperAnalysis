% Code below creates population summary for Significnat Saccade Direction Cells
% Written by Seth Konig
% Code does the following
% 1) Summarizes MRLs (mean resultant vector length) for place and non-place cells
% 2) Summarize circular non-uniformity p-value (biased by fixation count and firing rate)
% 3) Copies relevant figures for place cells to summary directory

%Code rechecked by SDK on 1/5/2017 & then re-written and rechecked again 5/3/2017

clar %clear,clc

%where to store spatial analysis figure copies
summary_directory =  'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PopulationFigures\Saccade Direction\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
saccades_with_spikes = [];
smval_deg = 6;%18 %9 degrees std

%---Values All Fixations out2out and in2in---%
all_mrls = []; %all observed MRLs (mean resultant vector length) ignoring out2in and in2out
all_mrl_pctiles = []; %observed MRLs shuffled percentile ignoring out2in and in2out

%---Values All Fixations out2out only---%
all_mlrs_out = []; %observed MRLs for out2out fixations only
all_mrl_out_pctiles = []; %observed MRLs shuffled percentile ignoring out2in and in2out

%---Other values---%
spatialness = []; %1 for place cell, 0 for non place cell
all_unit_names = {};
all_monkeys = []; %1s and 2s for monkeys
direction_cell_AP_location = []; %AP location of recorded place cell
place_cell_curves_in_min_out = [];


prefered_firing_rate_curves = [];
anti_prefered_firing_rate_curves = [];
ratio_prefered = [];
all_windows  = [];
normalized_prefered = [];
median_prefered_directions = [];

figure_dir = {};
all_dirs = []; %saccade directions for significant units
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML[
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'valid_trials',...
            'stability_attribute','cfg','hdr','data');
        
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
        
        disp(task_file(1:8))
        if exist([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat'])
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info')
            load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        twinad1 = 200;%presaccade window
        twinad2 = 400;%pos saccade window
        start_window = twinad1-min_fix_dur;%how early before saccade can you look for direction tuning
        end_window = twinad1+min_fix_dur+44;%how late after saccade can you look for direction tuning
        %minimum fixation duration + median saccade duration of 44 ms, some neurons have
        
        num_units = size(unit_stats,2);
        for unit =1:num_units
            if ~isnan(mrls.all_saccades(unit)) %if unit was processed
                
                %---Values All Fixations out2out and in2in---%
                all_mrls = [all_mrls mrls.all_saccades(unit)]; %all observed MRLs (mean resultant vector length) ignoring out2in and in2out
                all_mrl_pctiles = [all_mrl_pctiles mrls.all_saccades_shuffled_prctile(unit)];%observed MRLs shuffled percentile ignoring out2in and in2out
                
                all_dirs = [all_dirs saccade_direction{unit}];
                
                %---Values All Fixations out2out only---%
                all_mlrs_out = [all_mlrs_out mrls.out2out(unit)]; %observed MRLs for out2out fixations only
                all_mrl_out_pctiles = [all_mrl_out_pctiles mrls.out2out_shuffled_prctile(unit)]; %observed MRLs shuffled percentile ignoring out2in and in2out
                
                %---Other values---%
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monkey];
                if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1]; %place cell
                    
                    sac_firing = saccade_aligned_firing{unit};
                    in_out = sac_in_out{unit};
                    fix_ends = fixation_ends{unit};
                    
                    out2in_firing = sac_firing(in_out == 1,:);
                    o2i = nandens3(out2in_firing,smval,1000);
                    out2out_firing = sac_firing(in_out == 4,:);
                    o2o = nandens3(out2out_firing,smval,1000);
                    out2in_fixdurs = fix_ends(in_out == 1);
                    out2out_fixdurs = fix_ends(in_out == 4);
                    
                    [PKS,LOCS]= findpeaks(o2i-o2o,'MinPeakWidth',smval); %find peaks > 2*std in width
                    if isempty(PKS)
                        [PKS,LOCS]= findpeaks(o2i,'MinPeakWidth',smval); %find peaks > 2*std in width
                    end
                    if ~isempty(PKS)
                        LOCS = LOCS(PKS == max(PKS));
                        PKS = max(PKS);
                        median_fix_dur = median(fix_ends);
                        
                        if LOCS-twinad1 >  median_fix_dur
                            %can't really do much about it by can try
                            too_short = find(out2in_fixdurs < median_fix_dur);
                            out2in_firing(too_short,:) = [];
                            o2i = nandens3(out2in_firing,smval,1000);
                            
                            too_short = find(out2out_fixdurs < median_fix_dur);
                            out2out_firing(too_short,:) = [];
                            o2o = nandens3(out2out_firing,smval,1000);
                            
                            od =  o2i-o2o;
                            od = od-mean(od(1:twinad1));
                            od = od/max(od);
                            place_cell_curves_in_min_out = [place_cell_curves_in_min_out; od];
                        else%remove eye movements that are shorter than the peak
                            min_dur = LOCS-twinad1+smval; %add duration of smoothing param
                            too_short = find(out2in_fixdurs < min_dur);
                            out2in_firing(too_short,:) = [];
                            o2i = nandens3(out2in_firing,smval,1000);
                            
                            too_short = find(out2out_fixdurs < min_dur);
                            out2out_firing(too_short,:) = [];
                            o2o = nandens3(out2out_firing,smval,1000);
                            
                            od =  o2i-o2o;
                            od = od-mean(od(1:twinad1));
                            od = od/max(od);
                            place_cell_curves_in_min_out = [place_cell_curves_in_min_out; od];
                            
                        end
                    end
                else
                    spatialness = [spatialness 0]; %non place cell
                end
                direction_cell_AP_location = [direction_cell_AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                
                if  mrls.all_saccades_shuffled_prctile(unit) > 95
                    
                    window = all_direction_windows{unit};
                    all_windows = [all_windows {window}];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%---Grab Saccade Aligned Activity--%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    sac_aligned = saccade_aligned_firing{unit}; %fixation aligned firing rate
                    sac_dirs = saccade_directions{unit}; %saccade directions organized the same way as sac_algined
                    fix_starts = fixation_starts{unit};
                    fix_ends = fixation_ends{unit};
                    
                    %---Remove Counfounding Eye Movements for Spatially Modulated Neurons---%
                    %take eye data from fixations in2in or out2out since in2out or out2in
                    %could be biased by field location creating artificial direction tuning
                    % but only do this if spatial (both criterion)
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        sac_aligned = sac_aligned(sac_in_out{unit} == 2 | sac_in_out{unit} == 4,:); %fixations in2in or out2out
                        sac_dirs = sac_dirs(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for fixations in2in or out2out only
                        fix_starts = fix_starts((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
                        fix_ends = fix_ends((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
                    end
                    
                    %---Remove Fixations that are too short in duration---%
                    window_width = length(all_direction_windows{unit});
                    window_end = all_direction_windows{unit}(end);
                    window_start = all_direction_windows{unit}(1);
                    if window_end > end_window
                        window_end = all_direction_windows{unit}(end)-twinad1;
                        %remove fixations shorter than window as these could be
                        %contaminated by the next saccade
                        fixations_too_short = find(fix_ends < window_end);
                        sac_dirs(fixations_too_short) = [];
                        sac_aligned(fixations_too_short,:) = [];
                    end
                    if window_start < twinad1
                        fixations_too_short = find((twinad1+fix_starts) > window_start);
                        sac_dirs(fixations_too_short) = [];
                        sac_aligned(fixations_too_short,:) = [];
                    end
                    
                    [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,sac_aligned,window,sac_dirs);
                    binned_firing_rate_curves{1,unit} = mean_binned_firing_rate; %binned firing rates
                    [estimated_prefered_direction,prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
                        select_prefred_indeces(binned_firing_rate_curves{1,unit},degrees,sac_dirs,smval_deg,bin_deg);
                    
                    median_prefered_directions = [median_prefered_directions estimated_prefered_direction];
                    if strcmpi(task_file(1:8),'TO160515')
                        disp('now')
                    end
                    antipref = nandens3(sac_aligned(anti_prefered_dirs,:),smval,1000);
                    pref = nandens3(sac_aligned(prefered_dirs,:),smval,1000);
                    a = antipref;
                    a(a < 0.1) = mean(a);
                    a(a == 0) = 1;
                    ratio_prefered = [ratio_prefered; pref./a];
                    
                    %                     if max(pref./a) > 10
                    %                        disp('now')
                    %                     end
                    
                    if any(isnan(ratio_prefered(end,:)))
                        disp('now')
                    elseif sum(isinf(ratio_prefered(:))) > 0
                        disp('now')
                    end
                    
                    norm = pref-antipref;
                    norm = norm-mean(norm);
                    norm = norm/max(norm);
                    normalized_prefered = [normalized_prefered; norm];
                    
                    pref = pref-mean(pref);
                    pref = pref/max(abs(pref));
                    antipref = antipref-mean(antipref);
                    antipref = antipref/max(abs(antipref));
                    
                    prefered_firing_rate_curves = [prefered_firing_rate_curves; pref];
                    anti_prefered_firing_rate_curves = [anti_prefered_firing_rate_curves; antipref];
                else
                    sac_aligned = saccade_aligned_firing{unit}; %fixation aligned firing rate
                    sac_dirs = saccade_directions{unit}; %saccade directions organized the same way as sac_algined
                    sac_amps = saccade_amplitudes{unit}/24; %saccade directions organized the same way as sac_algined
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %---Determine if Unit is Significantly Modulated by Saccade Direction--%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %take eye data from fixations in2in or out2out since in2out or out2in
                    %could be biased by field location creating artificial direction tuning
                    % but only do this if spatial (both criterion)
                    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                        sac_aligned = sac_aligned(sac_in_out{unit} == 2 | sac_in_out{unit} == 4,:); %fixations in2in or out2out
                        sac_dirs = sac_dirs(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for fixations in2in or out2out only
                        sac_amps = sac_amps(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for fixations in2in or out2out only
                    end
                    %remove saccades that are too big since wont have many anyway
                    sac_aligned_amp = sac_aligned;
                    too_large = find(sac_amps > max_amplitude);
                    sac_aligned_amp(too_large,:) = [];
                    sac_amps(too_large)=[];
                    
                    window = all_direction_windows{unit};
                    spike_count = sum(sac_aligned_amp(:,window),2);
                    saccades_with_spikes = [saccades_with_spikes sum(spike_count > 0)];
                    
                end
            end
        end
    end
end
%%
clc
disp([num2str(nansum(spatialness)) ' place cells'])
disp([num2str(sum(all_mrl_pctiles > 95)) ' directionally modulated cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95)) ' directionally modulated cells for OUT2OUT fixations'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & all_mrl_pctiles > 95)) ' directionally modulated cells for OUT2OUT fixations & "all fixations"'])
disp('--------------------------------------------------------------')
disp([num2str(sum(all_mrl_pctiles > 95 & spatialness == 1)) ' directionally modulated place cells for "all fixations"'])
disp([num2str(sum(all_mrl_out_pctiles > 95 & spatialness == 1)) ' directionally modulated place cells for OUT2OUT fixations'])
disp('--------------------------------------------------------------')
disp(['Analyzed ' num2str(sum(spatialness == 1 & ~isnan(all_mrl_pctiles))) ' place cells'])
disp(['Analyzed ' num2str(sum(spatialness == 0 & ~isnan(all_mrl_pctiles))) ' non-place cells'])
disp(['Analyzed ' num2str(sum(~isnan(all_mrl_pctiles))) ' total cells'])
%%
%---Copy Relevant Figures to Summary Directory---%
for unit = 1:length(all_unit_names)
    if all_mrl_pctiles(unit) > 95
        sub_dir1 = '\Saccade Direction and Amplitude\';
        name1 = [all_unit_names{unit} '_Saccade_Direction_Analysis.png'];
        if spatialness(unit) == 1 %place cell
            if ~exist([summary_directory 'Place\'],'dir')
                mkdir([summary_directory 'Place\']);
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Place\' name1])
        elseif spatialness(unit) == 0 %non place cell
            if ~exist([summary_directory 'Non Place\'])
               mkdir([summary_directory 'Non Place\']); 
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Non Place\' name1])
        end
    end
end

%% Distribution of Saccade Directions
smval_deg =18; %9 degrees std

binned_dirs = zeros(1,360);
degrees = [-180:180];
for b = 1:length(degrees)-1
    binned_dirs(b) = sum(all_dirs < degrees(b+1) & all_dirs >= degrees(b));
end
degrees = degrees*pi/180;

bf = [binned_dirs(end-(3*smval_deg):end) binned_dirs binned_dirs(1:3*smval_deg)];%so don't get edege artifacts
binned_dirs = nandens(bf,smval_deg,'gauss',1); %smooth to display firing rate by saccade direction
binned_dirs = binned_dirs(3*smval_deg+2:end-3*smval_deg);%remove buffers
polarplot(degrees,[binned_dirs(end) binned_dirs])

%%
tm = -twin1:twin2-1;
figure
subplot(2,2,1)
hold on
plot(tm,nanmean(anti_prefered_firing_rate_curves),'r')
plot(tm,nanmean(prefered_firing_rate_curves),'g')
plot([tm(1) tm(end)],[0 0],'k')
plot([0 0],[-1 1],'k--')
hold off
xlabel('Time From Saccade Start (ms)')
ylabel('Normalized Firing Rate')
title('Population Average Mormalzed Firing Rates')
legend('Anti-prefered','Prefered')

subplot(2,2,2)
[~,mxi] = max(prefered_firing_rate_curves');
[~,si] = sort(mxi);
imagesc(tm,1:size(prefered_firing_rate_curves,1),prefered_firing_rate_curves(si,:))
xlabel('Time From Saccade Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Time Cell Plot')

subplot(2,2,3)
plot(tm,nandens3(ratio_prefered,20,1))
xlabel('Time From Saccade Start (ms)')
ylabel('Ratio')
title('Average Ratio of Prefered/Anti-Prefered')
box off


window_conc = zeros(1,twin1+twin2);
for w = 1:length(all_windows);
    window_conc(all_windows{w}) = window_conc(all_windows{w})+1;
end

subplot(2,2,4)
plot(tm,100*window_conc/length(all_windows))
xlabel('Time From Saccade Start (ms)')
ylabel('% of Units')
title('Distribution of 100 ms Window with Greatest Direction Modulation')
box off
%%
allvals = normalized_prefered;
allvals = allvals(allvals < 0);
stdvals = std(allvals);

%%
median_window = [];
window_len = [];
for n = 1:size(normalized_prefered,1)
    median_window(n) = median(all_windows{n});
    window_len(n) = length(all_windows{n});
    min_window = min(all_windows{n});
    max_window = max(all_windows{n});
    if min_window > 50
        min_window = min_window-50;
    else
        min_window = 1;
    end
    if max_window > twin2+twin1-50
        max_window = twin2+twin1;
    else
        max_window = max_window+50;
    end
    allvals(n,min_window:max_window) = NaN;
end

%%
%%
figure
subplot(1,2,1)
[~,mxi] = max(normalized_prefered');
[~,si] = sort(mxi);
imagesc(tm,1:size(normalized_prefered,1),normalized_prefered(si,:))
xlabel('Time From Saccade Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Time Cell Plot-sorted by max')
caxis([-stdvals 1])
axis square

sig_names = all_unit_names(all_mrl_pctiles > 95);
sig_names = sig_names(si);

subplot(1,2,2)
[~,si2] = sort(median_window);
imagesc(tm,1:size(normalized_prefered,1),normalized_prefered(si2,:))
xlabel('Time From Saccade Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Time Cell Plot-sorted by window')
caxis([-stdvals 1])
axis square

figure

spatial = normalized_prefered(spatialness(all_mrl_pctiles > 95) == 1,:);
subplot(2,2,1)
[~,mxi] = max(spatial');
[~,si] = sort(mxi);
imagesc(tm,1:size(spatial,1),spatial(si,:))
xlabel('Time From Saccade Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Spatial Neurons')
caxis([-stdvals 1])
axis square

window_conc = zeros(1,twin1+twin2);
spatial_windows = all_windows(spatialness(all_mrl_pctiles > 95) == 1);
for w = 1:length(spatial_windows);
    window_conc(spatial_windows{w}) = window_conc(spatial_windows{w})+1;
end

subplot(2,2,3)
plot(tm,100*window_conc/length(all_windows))
axis square
box off


nonspatial = normalized_prefered(spatialness(all_mrl_pctiles > 95) == 0,:);

subplot(2,2,2)
[~,mxi] = max(nonspatial');
[~,si] = sort(mxi);
imagesc(tm,1:size(nonspatial,1),nonspatial(si,:))
xlabel('Time From Saccade Start (ms)')
ylabel('Neuron #')
colormap('jet')
title('Non-Spatial Neurons')
caxis([-stdvals 1])
axis square

window_conc = zeros(1,twin1+twin2);
nonspatial_windows = all_windows(spatialness(all_mrl_pctiles > 95) == 0);
for w = 1:length(nonspatial_windows);
    window_conc(nonspatial_windows{w}) = window_conc(nonspatial_windows{w})+1;
end

subplot(2,2,4)
plot(tm,100*window_conc/length(all_windows))
axis square
box off
%%
pre_saccadic = [];
peri_saccadic = [];
post_saccadic = [];

for w = 1: length(all_windows)
    ind = all_windows{w};
    if all(ind < twinad1)
        pre_saccadic = [pre_saccadic w];
    elseif all(ind > twinad1+44)
        post_saccadic = [post_saccadic w];
    else
        peri_saccadic = [peri_saccadic w];
    end
end
%%
any_pre_sac = [];
for w = 1: length(all_windows)
    ind = all_windows{w};
    if sum(ind < twinad1) > smval/2 %smoothing window
        any_pre_sac = [any_pre_sac w];
    end
end
%%

figure
hold on
for m = 1:length(median_prefered_directions)
    polar([median_prefered_directions(m) median_prefered_directions(m)],[0 1],'k')
end
axis square
box off

figure
polarplot(0,0)
box off
%%
figure
hist(window_len,20)
hold on
plot([median(window_len) median(window_len)],[0 8],'r--')
hold off
xlabel('FWHM (ms)')
ylabel('Direction Cell Count')
box off
%%
sig_mrls = all_mrls(all_mrl_pctiles > 95);

figure
hist(sig_mrls,25)
hold on
plot([median(sig_mrls) median(sig_mrls)],[0 5],'r--')
hold off
xlabel('MRL')
ylabel('Saccade Direction Cell Count')
box off
title(['Median MRL: ' num2str(median(sig_mrls),2)])
%%
[~,mxi] = max(normalized_prefered');
[~,si] = sort(mxi);
sig_APs = direction_cell_AP_location(all_mrl_pctiles > 95);
figure
plot(sig_APs,mxi,'.k')

%%
vals = place_cell_curves_in_min_out(:,1:twinad1);
figure
[~,mxi] = max(place_cell_curves_in_min_out');
[~,place_order] = sort(mxi); %sort order by peak firing time
imagesc([-twinad1:twinad2-1],[1:size(place_cell_curves_in_min_out,1)],place_cell_curves_in_min_out(place_order,:))
hold on
plot([0 0],[0 size(place_cell_curves_in_min_out,1)],['w--']);
hold off
colormap('jet')
xlabel('Time From Saccade Start (ms)')
ylabel('View Cell #')
caxis([-std(vals(:)) 1])
colorbar
