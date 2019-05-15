function  Place_Cell_Fixation_Analysis(data_dir,figure_dir,session_data,predict_rt)
%written by Seth Konig 4/25/16 updated August 30, 2016
%Code determines/verifies that firing rate within the place field is great
%than firing rate outside of the place field. A common theme off
%view cells (based on an initial pass and a logical assumption one would make)
%is that they show strong perisaccadic modulation usually firing more after
%a fixation has started (hence based on viewing location). This code also
%analyzes firing rates during sequence trials for these neurons to
%determine if neuros show similar spatial representations during both
%tasks.
%
% Code rechecked bugs January 3, 2017 SDK
% Code re-rechecked for bugs February 3, 2017 SDK after significant
% modifications mostly around Maris method of mulitple comarpisons corrections


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Set important task parameters---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
place_cell_dir = [figure_dir 'Place Cells Fixation Analysis\']; %mother dir
best_place_cell_dir = [place_cell_dir 'Best Place Cells Fixation Analysis\']; %super conservative: skagg+corr 1/2
skagg_cell_dir99 = [place_cell_dir 'Skagg 99\'];%want to see what significant skagg scores > 99%-tile
skagg_cell_dir = [place_cell_dir 'Skagg only\'];%want to see what significant skagg scores shows
spatial_corr_dir = [place_cell_dir 'Spatial Correlation Only\']; %want to see what spatially consistent score shows
not_place_cell_dir = [place_cell_dir 'Not Spatial\'];%for neurons that don't pass any criterior
task = 'ListSQ';
colors ='rgbm';
shapes = 'xo';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task_file = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
try
    load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
catch
    disp('Spatial Analysis File not found continueing')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat']);

%Save as new variableso can reload later...kind of dumb but that was how it was written
absolute_fixationstats = fixationstats;
absolute_cfg = cfg;

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0; %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

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
[which_img,~,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%filter parameters
H = define_spatial_filter(filter_width);

%---Pre-allocate space for List Fixations Inside vs Outside Field Analysis---%
list_fixation_locked_firing = cell(1,num_units); %firing rate locked to fixations
saccade_direction = cell(1,num_units); %saccade direction organized the same as list_fixation_locked_firing
saccade_amplitude = cell(1,num_units); %saccade direction organized the same as list_fixation_locked_firing
area = NaN(1,num_units); %area of place field
all_place_field_matrix = cell(1,num_units); %place field location 1 for in 0 for out, NaN for no coverage
list_sig_times = cell(2,num_units); %%time points > 95% confidence interval, corrected for multiple comparisons using Maris method
in_out = cell(1,num_units); %organized same as list_fixation_locked_firing for fixations ....
%1) first fixation in: out-> in
%2) fixation in field but not first: in -> in
%3) first fixation out of field: in -> out
%4) fixation out of field but not first: out-> out

%---Pre-allocate space for Fixations during Sequence Trials Analysis---%
sequence_fixation_locked_firing = cell(4,num_units);%firing rate locked to fixations organized by item
sequence_sig_times = cell(2,num_units);  %shuffled 95% confidence interval for fixation in field vs out of field, row 1 97.5% row 2 2.5%
sequence_sig_times = cell(1,num_units);%time points > 95% confidence interval
all_which_sequence = cell(1,num_units); %which sequence
trial_nums = cell(1,num_units); %trial numbers
in_out_sequence = cell(1,num_units); %1 for sequence items inside 0 for suquence items outside NaN for items on border

stats_across_tasks = NaN(4,num_units);
%1) list peak location
%2) list peak
%3) sequence peak location
%4) sequence peak

for unit = 1:num_units

    if isempty(eyepos{unit})
        continue %no data for this neuron
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine List Firing rate Locked to Fixations Inside vs Outside Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Map---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
    place_field_matrix = define_place_field(firing_rate_map,imageX,imageY); %place field location
    all_place_field_matrix{unit} = place_field_matrix; %save place field location
    area(unit) = 100*nansum(nansum(place_field_matrix))/(imageX*imageY); %calculate area as % of screen size
    %note ignores incomplete coverage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Fixations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fix_in_out = NaN(1,3000); %firing rate for fixations inside or outside  of place field
    %1) first fixation in: out-> in
    %2) fixation in field but not first: in -> in
    %3) first fixation out of field: in -> out
    %4) fixation out of field but not first: out-> out
    fix_locked_firing = NaN(3000,(twin1+twin2)); %spike trains locked to fixations
    saccade_direction{unit} = NaN(1,3000);%saccade directions
    saccade_amplitude{unit} = NaN(1,3000);%saccade amplitudes
    
    fix_ind = 1; %fixation # so can track in variables above
    fixationstats = absolute_fixationstats; %reload because written over below
    cfg = absolute_cfg; %reload because written over below
    num_trials = length(cfg.trl);%number of trials
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
                    prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                    if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                        continue;
                    end
                    sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                    if sacamp < min_sac_amp %prior saccade is too small so ignore
                        continue
                    end
                    saccade_amplitude{unit}(fix_ind) = sacamp;

                    
                    prior_fix_in_out = [];
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
                    
                    %get firing rate locked to fixation
                    fixt = fixationtimes(1,f);%start of fixation
                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                    temp = zeros(1,twin1+twin2);
                    temp(fix_spikes) = 1;
                    fix_locked_firing(fix_ind,:) = temp;
                    
                    %reflix y otherwise everything is flipped
                    fixy = fixations(2,f);
                    last_fixy = fixations(2,f-1);
                    saccade_direction{unit}(fix_ind) = atan2d(fixy-last_fixy,fixx-last_fixx);

                    fix_ind = fix_ind+1;
                end
            end
        end
    end
    
    %remove excess NaNs
    saccade_direction = laundry(saccade_direction);
    saccade_amplitude = laundry(saccade_amplitude);
    fix_locked_firing = laundry(fix_locked_firing);
    fix_in_out = laundry(fix_in_out);
    
    %store variables across units for later access
    list_fixation_locked_firing{unit} = fix_locked_firing;
    in_out{unit} = fix_in_out;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%---Determine if Unit Passes 95% for both Skaggs and Stability--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(spatial_info.shuffled_rate_prctile(unit))
        continue %spatial analysis wasn't run on this unit
    elseif (spatial_info.shuffled_rate_prctile(unit) < 95) && (spatial_info.spatialstability_halves_prctile(unit) < 95)
        continue %unit not spatial
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neuron Fires Faster for Fixations in the Field vs Out of Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %---Permutation/Resampling Analysis for Significance---%
    %for fixations out->in vs out->out
    firing_rate_in_out = fix_locked_firing(fix_in_out == 1 | fix_in_out == 4,:);
    [~,firing_rate_curve] = nandens(firing_rate_in_out,smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
    in_or_out = fix_in_out(fix_in_out == 1 | fix_in_out == 4);
    all_curves = NaN(numshuffs,twin1+twin2);
    
    %for fixations ?->in vs ?->out in case not enough samples for the one above
    in_or_out2 = fix_in_out;
    [~,all_firing_rate_curve] = nandens(fix_locked_firing,smval,'gauss',Fs,'nanflt'); %trial-by-trial smoothed
    all_curves2 = NaN(numshuffs,twin1+twin2);
    
    parfor shuff = 1:numshuffs;
        %for fixations out->in vs out->out
        ind = randperm(length(in_or_out));
        shuff_in_or_out = in_or_out(ind);%randomly assign fixations as in or out
        shuff_in_curve = nanmean(firing_rate_curve(shuff_in_or_out == 1,:)); %out->in
        shuff_out_curve = nanmean(firing_rate_curve(shuff_in_or_out == 4,:));%out->out
        all_curves(shuff,:) = shuff_in_curve-shuff_out_curve; %calculate shuffled difference
        
        %for fixations ?->in vs ?->out in case not enough samples for the one above
        ind2 = randperm(length(in_or_out2));
        shuff_in_or_out = in_or_out2(ind2);
        shuff_in_curve = nanmean(all_firing_rate_curve(shuff_in_or_out == 1 | shuff_in_or_out == 2,:)); %?->in
        shuff_out_curve = nanmean(all_firing_rate_curve(shuff_in_or_out == 3 | shuff_in_or_out == 4,:));%?->out
        all_curves2(shuff,:) = shuff_in_curve-shuff_out_curve; %calculate shuffled difference
    end
    
    %for fixations out->in vs out->out
    in_curve = nanmean(firing_rate_curve(in_or_out == 1,:));  %firing rate curve for fixations in the field
    list_in_curve = in_curve; %firing rate curve for fixations in the field
    out_curve = nanmean(firing_rate_curve(in_or_out == 4,:)); %firing rate fo fixations out of the filed
    observed_diff = in_curve-out_curve; %observed difference in firing rate 
    [~,list_sig_times{1,unit}] = cluster_level_statistic(observed_diff,all_curves,1,smval); %multiple comparision corrected significant indeces
    
    %for fixations ?->in vs ?->out in case not enough samples for the one above    
    in_curve = nanmean(all_firing_rate_curve(in_or_out2 == 1 | in_or_out2 == 2,:));
    out_curve = nanmean(all_firing_rate_curve(in_or_out2 == 3 | in_or_out2 == 4,:));
    observed_diff = in_curve-out_curve; %observed difference in firing rate 
    [~,list_sig_times{2,unit}] = cluster_level_statistic(observed_diff,all_curves2,1,smval); %multiple comparision corrected significant indeces
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Peak and Time of Peak Firing Rate aligned to Fixation----%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    in_curve = list_in_curve; %firing rate curve for fixations in the field
    in_curve = in_curve-nanmean(in_curve(1:twin1)); %subtract out of field firing rate
    in_curve = in_curve/max(in_curve); %normalize to maximum firing rate
    [PKS,LOCS]= findpeaks(in_curve,'MinPeakWidth',smval); %find peaks > 2*std in width
    %find peaks that were in significant period and remove those that weren't
    rmv = [];
    for l = 1:length(LOCS);
        if list_sig_times{1,unit}(LOCS(l)) == 0 && list_sig_times{2,unit}(LOCS(l)) == 0  %is it at a reliable time
            rmv = [rmv l];
        end
    end
    LOCS(rmv) = [];
    PKS(rmv) = [];
    if ~isempty(LOCS)
        %remove peaks less than 2/3 the max
        LOCS(PKS < 0.66) = [];
    end
    if ~isempty(LOCS)
        LOCS = LOCS(1);%if multiple peaks take first peak
        stats_across_tasks(1,unit) = LOCS; %time of peak
        stats_across_tasks(2,unit) = list_in_curve(LOCS); %peak firing rate
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Make Plot for List Fixations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t = -twin1:twin2-1;
    maxfr = prctile(firing_rate_map(:),97.5);
    
    figure
    
    %---plot firing rate map for all images---%
    subplot(2,3,1)
    h = imagesc(firing_rate_map);
    set(h,'alphadata',~isnan(firing_rate_map));
    title('All images')
    axis off
    axis equal
    colorbar
    colormap('jet')
    clim = caxis;
    caxis([clim(1) maxfr])
    
    title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
    if spatial_info.shuffled_rate_prctile(unit) >= 95;
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
    
    
    %---draw place field---%
    subplot(2,3,4)
    imagesc(place_field_matrix)
    axis off
    axis equal
    title(sprintf(['Place Field Location \n Area = ' num2str(area(unit),2)]));
    
    %---Raster: Fixations in->out vs out->out---%
    subplot(2,3,2)
    [trial,time] = find(fix_locked_firing(fix_in_out == 4,:) == 1); %out2out
    plot(time-twin1,(trial),'.b')
    hold on
    if ~isempty(trial)
        b4 = max(trial);
    else
        b4 = 0;
    end
    [trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1); %out2in
    trial = trial+b4;
    plot(time-twin1,(trial),'.r')
    if ~isempty(trial)
        ylim([0 max(trial)])
    else
        ylim([0 b4])
    end
    box off
    ylabel('Occurence #')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    xlabel('Time from Fixation Start (ms)')
    title('Fixation Aligned Rasters')
    
    %---Firing Rate Curve: Fixations in->out vs out->out---%
    subplot(2,3,5)
    hold on
    dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
    dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
    if ~isnan(stats_across_tasks(1,unit))
        plot(stats_across_tasks(1,unit)-twin1,list_in_curve(stats_across_tasks(1,unit)),'*k')
    end
    yl = ylim;
    if yl(1) < 0
        yl(1) = 0;
        ylim(yl);
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    gaps = findgaps(find(list_sig_times{1,unit}));
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
    xlim([-twin1 twin2]);
    hold off
    xlabel('Time from Fixation Start (ms)')
    ylabel('Firing Rate (Hz)')
    legend('out->in','out->out','Location','NorthWest')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
        ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
    
    %---Raster: Fixations ?->in vs ?-> out---%
    subplot(2,3,3)
    [trial,time] = find(fix_locked_firing(fix_in_out == 3 | fix_in_out == 4,:) == 1); %?2out
    plot(time-twin1,(trial),'.b')
    hold on
    if ~isempty(trial)
        b4 = max(trial);
    else
        b4 = 0;
    end
    [trial,time] = find(fix_locked_firing(fix_in_out == 1 | fix_in_out == 2,:) == 1);%?2in
    trial = trial+b4;
    plot(time-twin1,(trial),'.r')
    if ~isempty(trial)
        ylim([0 max(trial)])
    else
        ylim([0 b4])
    end
    ylabel('Occurence #')
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    xlabel('Time from Fixation Start (ms)')
    title('Fixation Aligned Rasters')
    xlim([-twin1 twin2])
    hold off
    
    %---Firing Rate Curve: Fixations ?->in vs ?-> out---%
    subplot(2,3,6)
    hold on
    dofill(t,fix_locked_firing(fix_in_out == 1 | fix_in_out == 2,:),'red',1,smval);%?->in
    dofill(t,fix_locked_firing(fix_in_out == 3 | fix_in_out == 4,:),'blue',1,smval);%?->out
    legend('?->in','?->out','Location','NorthWest')
    yl = ylim;
    if yl(1) < 0
        yl(1) = 0;
        ylim(yl);
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    gaps = findgaps(find(list_sig_times{2,unit}));
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
    box off
    hold off
    xlabel('Time from Fixation Start (ms)')
    ylabel('Firing Rate (Hz)')
    xlim([-twin1 twin2])
    set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
    title(sprintf(['n_{?->in} = ' num2str(sum(fix_in_out == 1 | fix_in_out == 2))...
        ', n_{?->out} = ' num2str(sum(fix_in_out == 3 | fix_in_out == 4))]));
    
    num_trials_str = [' n_{img} = ' num2str(size(spike_times{unit},1))];
    if multiunit(unit)
        multi_str = 'Multiunit ';
    else
        multi_str = ' ';
    end
    subtitle(['Spatial Plots' num_trials_str multi_str ' ' task_file(1:8) ' ' unit_stats{1,unit}]);
    
    if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)% ...
        save_and_close_fig(best_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 95) && ((spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
            || (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95));
        save_and_close_fig(spatial_corr_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 99)
        save_and_close_fig( skagg_cell_dir99,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    elseif (spatial_info.shuffled_rate_prctile(unit) > 95)
        save_and_close_fig(skagg_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    else
        save_and_close_fig(not_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_fixation_analysis']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Spatial Representations are Similar in Sequence Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %preallocate space and parallel structure of cfg
    successful_sequence_trials = NaN(1,length(cfg.trl)); %sequence trial number
    which_sequence = NaN(1,length(cfg.trl)); %which sequence 1 or 2
    for t = 1:length(cfg.trl);
        if sum(cfg.trl(t).allval == 3) == 6; %in which sequence trials were rewarded
            which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
            successful_sequence_trials(t) = t;
        end
    end
    successful_sequence_trials = laundry(successful_sequence_trials); %remove nans
    which_sequence = laundry(which_sequence); %remove nans
    all_which_sequence{unit} = which_sequence; %store for later
    num_trials = length(successful_sequence_trials); %number of sequence trials
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data locked to trial events---%%%
    fixationstats = absolute_fixationstats; %reload because will be written over below
    cfg = absolute_cfg; %reload because will be written over below
    fixationstats = fixationstats(successful_sequence_trials);
    cfg.trl = cfg.trl(successful_sequence_trials);
    
    saccade_start_time = NaN(num_trials,4);%when did saccade to item start
    fixation_start_time = NaN(num_trials,4);%when did fixation on item start
    time_to_fixation = NaN(num_trials,4); %how fast did they get to it the first time around
    fixation_duration = NaN(num_trials,4); %fixation duration
    for trial = 1:num_trials
        locs = sequence_locations{which_sequence(trial)};
        
        %Must convert to DVA for this analysis
        locs(1,:) = (locs(1,:)-imageX/2)/24; %convert to dva
        locs(2,:) = (locs(2,:)-imageY/2)/24; %convert to dva
        fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-imageX/2)/24;
        fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-imageY/2)/24;
        
        trial_codes = cfg.trl(trial).allval;
        trial_codes(trial_codes == 100)= 0; %have to set to 0 because code below runs behavior only data too
        trial_codes(1) = 100;%eye data starts for recording right away
        event_times = cfg.trl(trial).alltim;
        
        trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
            trial_codes,event_times,true); %runs all eye data analysis
        
        time_to_fixation(trial,:) = trialdata.t2f;
        fixation_duration(trial,:) = trialdata.fixation_duration;
        
        fixation_numbers = trialdata.fixationnums; %fixation number for each item
        fixationtimes = fixationstats{trial}.fixationtimes;
        saccadetimes = fixationstats{trial}.saccadetimes;
        for item = 1:4
            if ~isnan(fixation_numbers(item))
                fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
                saccadeind = find(saccadetimes(2,:)+1 ==  fixation_start_time(trial,item));
                if ~isempty(saccadeind)
                    saccade_start_time(trial,item) = saccadetimes(1,saccadeind);
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate Firing Rate Locked to Eye Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c = 1:4
        sequence_fixation_locked_firing{c,unit} = NaN(num_trials,twin1+twin2); %spike trains aligned to fixations
        trial_nums{c,unit} = NaN(1,num_trials); %trial number
    end
    
    for trial = 1:num_trials
        if successful_sequence_trials(trial) >= valid_trials(1,unit) && successful_sequence_trials(trial) <= valid_trials(2,unit)
            spikes = find(data(unit).values{successful_sequence_trials(trial)}); %trial spike trains
            for c = 1:4;
                %                 if time_to_fixation(trial,c) < predict_rt %remove predictive saccades for now
                %                     continue
                %                 end
                fixt = fixation_start_time(trial,c);
                sact = saccade_start_time(trial,c);
                
                %if no saccade detected on item want to check if previous
                %trial has saccade/begining of fixation
                if isnan(sact) && isnan(fixt) %both could not be located
                    continue
                elseif ~isnan(sact) && isnan(fixt) %could not locate fixation or saccade
                    error(['Could not identify fixation but could locate saccade! Trial#:' num2str(trial)])
                else
                    fix_spikes = spikes(spikes > fixt-twin1 & spikes <= fixt+twin2)-fixt+twin1;
                    temp = zeros(1,twin1+twin2);
                    temp(fix_spikes) = 1;
                    sequence_fixation_locked_firing{c,unit}(trial,:) = temp;
                    trial_nums{c,unit}(trial) = trial;
                end
            end
        end
    end
    sequence_fixation_locked_firing(:,unit) = laundry(sequence_fixation_locked_firing(:,unit));
    trial_nums(:,unit) = laundry(trial_nums(:,unit));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine Relationship Between List and Sequence Firing Rates---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Determine if any of the items are inside the matrix or not?
    sequence_inside = NaN(2,4);
    for c = 1:4
        for seq = 1:2
            yloc = imageY-sequence_locations{seq}(2,c);
            xloc = sequence_locations{seq}(1,c);
            
            %---Simple version...is item in or out of field---%
            %             if place_field_matrix(yloc,xloc) == 1 %item is in field
            %                 sequence_inside(seq,c) = 1;
            %             elseif place_field_matrix(yloc,xloc) == 0 %item is out of the field
            %                 sequence_inside(seq,c) = 0;
            %             else %coverage unknown
            %                 sequence_inside(seq,c) = NaN;
            %             end
            
            %---Less Simple Version...is item on border of field---%
            if place_field_matrix(yloc,xloc) == 1 %item is in field
                %then check if item is on border of field, if yes don't
                %count. Check if field extends at least 0.5 dva, should be
                %effected by coverage on edges items are at least 3.5 dva
                %away from border of image. Median eye tracking error (fixation accuracy) 
                %on this task ~0.56 dva so should be ok with 0.5 and most
                %fixations within 1 dva
                if place_field_matrix(yloc-12,xloc) == 1 && place_field_matrix(yloc-12,xloc-12) == 1 &&...
                        place_field_matrix(yloc-12,xloc+12) == 1 && place_field_matrix(yloc+12,xloc) == 1 && ...
                        place_field_matrix(yloc+12,xloc-12) == 1 && place_field_matrix(yloc+12,xloc+12) == 1 && ...
                        place_field_matrix(yloc,xloc+12) == 1 && place_field_matrix(yloc,xloc-12) == 1
                    sequence_inside(seq,c) =1;
                else
                    sequence_inside(seq,c) = NaN; %don't want to use border for any category
                end
            else %check if item outside is also close to the border
                if place_field_matrix(yloc-12,xloc) == 1 || place_field_matrix(yloc-12,xloc-12) == 1 ||...
                        place_field_matrix(yloc-12,xloc+12) == 1 || place_field_matrix(yloc+12,xloc) == 1 || ...
                        place_field_matrix(yloc+12,xloc-12) == 1 || place_field_matrix(yloc+12,xloc+12) == 1 || ...
                        place_field_matrix(yloc,xloc+12) == 1 || place_field_matrix(yloc,xloc-12) == 1
                    sequence_inside(seq,c) =NaN; %don't want to use border for any category
                else
                    sequence_inside(seq,c) = 0;
                end
            end
        end
    end
    
    if any(sequence_inside(:) == 1) && any(sequence_inside(:) == 0)
        seq_in_out = [];
        fixation_firing = [];
        for c = 1:4
            for seq = 1:2
                fixation_firing = [fixation_firing; ...
                    sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq,:)];
                if  sequence_inside(seq,c) == 1;
                    seq_in_out = [ seq_in_out ones(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
                elseif isnan(sequence_inside(seq,c))
                    seq_in_out = [ seq_in_out NaN(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
                else
                    seq_in_out = [ seq_in_out zeros(1,sum(which_sequence(trial_nums{c,unit}) == seq))];
                end
            end
        end
        in_out_sequence{unit} = seq_in_out;
        [~,firing_rate_curves] = nandens(fixation_firing,smval,'gauss',Fs,'nanflt'); %trial by trial smoothing
        all_curves = NaN(numshuffs,twin1+twin2);
        parfor shuff = 1:numshuffs;
            ind = randperm(length(seq_in_out));
            shuff_in_or_out = seq_in_out(ind);
            shuff_in_curve = nanmean(firing_rate_curves(shuff_in_or_out == 1,:));
            shuff_out_curve = nanmean(firing_rate_curves(shuff_in_or_out == 0,:));
            all_curves(shuff,:) = shuff_in_curve-shuff_out_curve;
        end
        in_curve = nandens(fixation_firing(seq_in_out == 1,:),smval,'gauss',Fs,'nanflt');
        out_curve = nandens(fixation_firing(seq_in_out == 0,:),smval,'gauss',Fs,'nanflt');
        observed_diff = in_curve-out_curve; %observed difference in firing rate 
        [~,sequence_sig_times{unit}] = cluster_level_statistic(observed_diff,all_curves2,2,smval); %multiple comparision corrected significant indeces
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Determine Peak and Time of Peak Firing Rate aligned to Fixation----%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        seq_in_curve = in_curve; %firing rate curve for fixations in the field
        in_curve = in_curve-nanmean(in_curve(1:twin1)); %subtract out of field firing rate
        in_curve = in_curve/max(in_curve); %normalize to maximum firing rate
        [PKS,LOCS]= findpeaks(in_curve,'MinPeakWidth',smval); %find peaks > 2*std in width
        %find peaks that were in significant period and remove those that weren't
        rmv = [];
        for l = 1:length(LOCS);
            if sequence_sig_times{unit}(LOCS(l)) == 0 %is it at a reliable time?
                rmv = [rmv l];
            end
        end
        LOCS(rmv) = [];
        PKS(rmv) = [];
        if ~isempty(LOCS)
            %remove peaks less than 2/3 the max
            LOCS(PKS < 0.66) = [];
        end
        if ~isempty(LOCS)
            LOCS = LOCS(1);
            stats_across_tasks(3,unit) = LOCS; %time of peak
            stats_across_tasks(4,unit) = seq_in_curve(LOCS); %peak firing rate
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Plot Sequence Firing Rates for Items Inside vs Outside Field---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t = -twin1:twin2-1;
        
        figure
        
        %---Plot firing rate map---%
        a = subplot(2,3,1);
        h = imagesc(firing_rate_map);
        set(h,'alphadata',~isnan(firing_rate_map));
        axis off
        axis equal
        colormap(a,'jet')
        colorbar
        caxis([clim(1) maxfr])
        title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
        if spatial_info.shuffled_rate_prctile(unit) >= 95;
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
        
        %---Plot Place Field---%
        b = subplot(2,3,2);
        imagesc(place_field_matrix);
        hold on
        for c = 1:4
            for seq = 1:2
                plot(sequence_locations{seq}(1,c),imageY-sequence_locations{seq}(2,c),[colors(c) shapes(seq)],'markersize',16)
            end
        end
        hold off
        xlim([0 800])
        ylim([0 240])
        axis equal
        axis off
        colormap(b,'gray')
        title(sprintf(['Place Field Location \n Area = ' num2str(area(unit),2) '%%']));
        
        %---Plot Firing Rate Curves for Each Item---%
        leg = {};
        subplot(2,3,3)
        hold on
        for seq = 1:2
            for c = 1:4
                [Dmean,~] = nandens (sequence_fixation_locked_firing{c,unit}(which_sequence(trial_nums{c,unit}) == seq ,:), smval, 'gauss',Fs);
                plot(t(1:10:end),Dmean(1:10:end),[colors(c) shapes(seq) '-']);%plot every 10 so can see visually
                leg = [leg {[num2str(c) shapes(seq)]}];
            end
        end
        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate (Hz)')
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlim([-twin1 twin2])
        title('Firing Rate Curves for Each Item')
        legend(leg,'Location','NorthEastOutside')
        
        %---Plot Firing Rate Curves for Suquence Trials---%
        subplot(2,3,4)
        hold on
        dofill(t,fixation_firing(seq_in_out == 1,:),'red',1,smval);
        dofill(t,fixation_firing(seq_in_out == 0,:),'blue',1,smval);
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        gaps = findgaps(find(sequence_sig_times{unit}));
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
        plot([0 0],[yl(1) yl(2)],'k--')
        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate (Hz)')
        legend('Items Inside','Items Outside','Location','NorthWest')
        if ~isnan(stats_across_tasks(3,unit))
            plot(stats_across_tasks(3,unit)-twin1,stats_across_tasks(4,unit),'*k')
            title(['Sequence Trials: peak of ' num2str(stats_across_tasks(4,unit),3) 'Hz @ ' ...
                num2str(stats_across_tasks(3,unit)-twin1) ' ms']);
        else
            title(['Sequence Trials: No peak, max firing of ' num2str(max(seq_in_curve),3) 'Hz']);
        end
        hold off
        
        %---Plot Firing Rate Curves for Image Trials---%
        subplot(2,3,5)
        hold on
        dofill(t(1:twin1+twin2),list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'red',1,smval);%out-> in
        dofill(t(1:twin1+twin2),list_fixation_locked_firing{unit}(in_out{unit} == 4,:),'blue',1,smval);%out->out
        plot(stats_across_tasks(1,unit)-twin1,stats_across_tasks(2,unit),'*k')
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        gaps = findgaps(find(list_sig_times{1,unit}));
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
        box off
        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate (Hz)')
        legend('out->in','in->in','Location','NorthWest')
        if ~isnan(stats_across_tasks(1,unit))
            plot(stats_across_tasks(1,unit)-twin1,stats_across_tasks(2,unit),'*k')
            title(['Image Trials: peak of ' num2str(stats_across_tasks(2,unit),3) 'Hz @ ' ...
                num2str(stats_across_tasks(1,unit)-twin1) ' ms']);
        else
            title(['Sequence Trials: No peak, max firing of ' num2str(max(list_in_curve),3) 'Hz']);
        end
        hold off
        
        %---Plot Firing Rate Curves for Suquence vs Image Trials---%
        subplot(2,3,6)
        hold on
        dofill(t(1:twin1+twin2),list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'black',1,smval);%list out-> in
        dofill(t,fixation_firing(seq_in_out == 1,:),'green',1,smval);%sequence
        yl = ylim;
        if yl(1) < 0
            yl(1) = 0;
            ylim(yl);
        end
        plot([0 0],[yl(1) yl(2)],'k--')
        hold off
        xlabel('Time from Fixation Start (ms)')
        ylabel('Firing Rate (Hz)')
        legend('Images','Sequences','Location','NorthWest')
        
        if ~isnan(stats_across_tasks(1,unit)) && ~isnan(stats_across_tasks(3,unit))
            title(['Contextual Gain: ' num2str(100*(stats_across_tasks(2,unit)/stats_across_tasks(4,unit)),3) '%'])
        end
        
        subtitle([task_file(1:8) ' ' unit_stats{1,unit}]);
        
        if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95) %...
            save_and_close_fig(best_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 95) && ((spatial_info.spatialstability_halves_prctile(1,unit) > 95) ...
                || (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95));
            save_and_close_fig(spatial_corr_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 99)
            save_and_close_fig( skagg_cell_dir99,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        elseif (spatial_info.shuffled_rate_prctile(unit) > 95)
            save_and_close_fig(skagg_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        else
            save_and_close_fig(not_place_cell_dir,[task_file(1:end-11) '_' unit_stats{1,unit} '_place_cell_Sequence_InSideOutside']);
        end
    end
end

save([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'],...
    'twin1','twin2','smval','list_fixation_locked_firing','in_out',...
    'task_file','area','stats_across_tasks','unit_stats','sequence_fixation_locked_firing',...
    'all_place_field_matrix','trial_nums','numshuffs','in_out_sequence',...
    'all_which_sequence','saccade_direction','list_sig_times','sequence_sig_times',...
    'saccade_amplitude')
end