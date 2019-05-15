% Population Visual Response
% Code below creates population summary for Significnat Visually Reponsive Neurons
% Written by Seth Konig  written Seth Konig 9/1/16, updated 1/16/2017
% Code does the following
% 1) Summarizes visual responses to images for short and long windows
% 2) Determines if neurons may be sequentially organized in time
% 3) Determines whether place cells are also visually responsive
% 4) Determines if visually responsive neurons are also modulated by novel/repeat images
% 5) Tracks AP location, unit counts, and which monkey (not currently used)
% 6) Copies relevant figures to summary directory

%Code rechecked by SDK on 1/16/2017

clar %clear,clc
%where to store spatial analysis figure copies
summary_directory = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PopulationFigures\Visual Response\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size

%---Misc. Parameters---%
monkey_all_unit_count = zeros(2,2);%row 1 place row 2 non-place, column by monkey
all_multi_unit_count = zeros(1,2); %all_multi_unit_countunits
all_unit_names = {}; %place cell unit names
all_monkeys = []; %1s and 2s
AP_location = []; %AP location of recorded place cell

%---Unit Test Significance---%
visual_response_stats_short = []; %short 1 second window, 1 for sig, 0 for not sig
visual_response_stats_long = []; %long 5 second window, 1 for sig, 0 for not sig
visual_response_stats_short_sliding = []; %sliding window analysis for short 1 second window, 1 for sig, 0 for not sig
visual_response_stats_long_sliding = []; %sliding window analysis for long 5 second window, 1 for sig, 0 for not sig
spatialness = []; %1 for place cells, 0 for non place cells
all_list_peak_times = [];%peak firing times for view cells
fixation_on_cross_status = []; %1 for sig, 0 for not sig
image_off_status = [];%1 for sig, 0 for not sig
memory_short = [];%significant indeces within the first 1 second window
memory_long = [];%significant indeces for 5 second window

%---Firing Rate Curves---%
cross_fixation_firing_rates = [];%fixation on across
image_on_firing_rates = []; %image on 1 second window
long_image_on_firing_rates = [];%image on 5 second window
image_off_firing_rates = []; %image off


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
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PW Recording Files\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PW Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = 'P:\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
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
        if exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file') %want to remove later
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']) %visual response analysis
            load([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results']) %memory visual response analysis
            load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat'],'spatial_info') %spatial analysis
            load([data_dir task_file(1:8) '-Place_Cell_Analysis.mat'])

        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        for unit = 1:num_units
            if ~isempty(time_lock_firing{unit,1})
                
                %---Misc. Parameters---%
                monkey_all_unit_count(1,monk) = monkey_all_unit_count(1,monk)+1; %unit count
                ll_multi_unit_count = all_multi_unit_count+multiunit(unit); %multiunit count
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}]; %unit name
                all_monkeys = [all_monkeys monk]; %1s and 2s
                AP_location = [AP_location chamber_zero(1)+ session_data{sess}.location(1)]; %AP location of recorded place cell
                
                
                %---Unit Test Significance and Firing Rate Curves---%
                if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1]; %place cell
                    all_list_peak_times = [all_list_peak_times stats_across_tasks(1,unit)];%time of peak firing rate of place cells for out-> in fixations

                else
                    spatialness = [spatialness 0]; %non place cell
                     all_list_peak_times = [all_list_peak_times NaN];
                end
                
                %---for fixation on cross hair---%
                if (epoch_data.rate_prctile(unit,2) > 95) && (epoch_data.temporalstability_prctile(unit,2) > 95)
                    fixation_on_cross_status = [fixation_on_cross_status 1]; %significant
                    firing_rate = time_lock_firing{unit,2}(:,1:twin1+twin3); %for now only look at first 500 ms after fixation
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    cross_fixation_firing_rates = [cross_fixation_firing_rates; firing_rate];
                else
                    fixation_on_cross_status = [fixation_on_cross_status 0]; %not significan
                end
                
                
                %---for image on short 1 second window---%
                if (epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95)
                    visual_response_stats_short = [visual_response_stats_short 1]; %signficant
                    firing_rate = time_lock_firing{unit,3};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_on_firing_rates = [image_on_firing_rates ; firing_rate];
                else
                    visual_response_stats_short = [visual_response_stats_short 0]; %not signficant
                end
                
                
                %---for long image on 5 second period---%
                if (epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)
                    visual_response_stats_long = [ visual_response_stats_long 1];%signficant
                    firing_rate = time_lock_firing{unit,5};
                    [firing_rate,~]= nandens(firing_rate,smval2,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    long_image_on_firing_rates = [long_image_on_firing_rates ; firing_rate];
                else
                    visual_response_stats_long = [ visual_response_stats_long 0];%not signficant
                end
                
                %---for image off---%
                if (epoch_data.rate_prctile(unit,4) > 95) && (epoch_data.temporalstability_prctile(unit,4) > 95)
                    image_off_status = [image_off_status 1];%signficant
                    firing_rate = time_lock_firing{unit,4};
                    [firing_rate,~]= nandens(firing_rate,smval,'gauss',Fs,'nanflt');%going to smooth slightl lesss than before
                    firing_rate = firing_rate-nanmean(firing_rate);
                    if max(firing_rate) > abs(min(firing_rate)) %normalize to max
                        firing_rate = firing_rate/max(abs(firing_rate));
                    else%normalize to min %could have some neurons that only show supression
                        firing_rate = firing_rate/min(firing_rate);
                    end
                    image_off_firing_rates = [image_off_firing_rates ; firing_rate];
                else
                    image_off_status = [image_off_status 0];%not signficant
                end
                
                %---Memory Responses---%
                if visual_response_stats_short(end) == 1 %then visual responsive within 1 second window
                    if sum(sig_short{unit}) > 0
                        gaps = findgaps(find(sig_short{unit}));
                        rmv = [];
                        for g = 1:size(gaps,1)
                            gp = gaps(g,:);
                            gp(gp == 0) = [];
                            if any(gp < twin1)
                                if all(gp < twin1)
                                    rmv = [rmv g];
                                else
                                    gp(gp < twin1) = [];
                                    if length(gp) < smval
                                        rmv = [rmv g];
                                    end
                                end
                            end
                        end
                        gaps(rmv,:) = [];
                        if sum(sum(gaps)) > 0
                            memory_short = [memory_short 1];
                        else
                            memory_short = [memory_short 0];
                        end
                    else
                        memory_short = [memory_short 0];
                    end
                else
                    memory_short = [memory_short NaN];
                end
                
                if visual_response_stats_long(end) == 1 %then visual responsive within 5 second window
                    if sum(sig_long{unit}) > 0
                        gaps = findgaps(find(sig_long{unit}));
                        rmv = [];
                        for g = 1:size(gaps,1)
                            gp = gaps(g,:);
                            gp(gp == 0) = [];
                            if any(gp < twin1)
                                if all(gp < twin1)
                                    rmv = [rmv g];
                                else
                                    gp(gp < twin1) = [];
                                    if length(gp) < smval2
                                        rmv = [rmv g];
                                    end
                                end
                            end
                        end
                        gaps(rmv,:) = [];
                        if sum(sum(gaps)) > 0
                            memory_long = [memory_long 1];
                        else
                            memory_long = [memory_long 0];
                        end
                    else
                        memory_long = [memory_long 0];
                    end
                else
                    memory_long = [memory_long NaN];
                end
                
                
                %---Sliding Window Analysis Simialr to Jutras and Buffalo---%
                if sum(sig_visual_response(unit,:)) == 0%nothing is significant
                    visual_response_stats_short_sliding = [visual_response_stats_short_sliding 0];
                    visual_response_stats_long_sliding = [visual_response_stats_long_sliding 0];
                else %something is signfiicant so check short and long times
                    if sum(sig_visual_response(unit,1:twin1+twin2)) > 0 %then short has significant response
                        visual_response_stats_short_sliding = [visual_response_stats_short_sliding 1];
                    else %not significant
                        visual_response_stats_short_sliding = [visual_response_stats_short_sliding 0];
                    end
                    if sum(sig_visual_response(unit,twin1+twin2:end)) > 0 %then long has significant response
                        visual_response_stats_long_sliding = [visual_response_stats_long_sliding 1];
                    else %then not significant
                        visual_response_stats_long_sliding = [visual_response_stats_long_sliding 0];
                    end
                end
            end
        end
    end
end
%% Display Summary Results
clc
disp([num2str(nansum(spatialness)) ' place cells'])
disp([num2str(nansum(spatialness == 1 & isnan(all_list_peak_times))) ' place cells have no peaks'])
disp([num2str(sum(visual_response_stats_short == 1)) ' Visual Response 1 second window'])
disp([num2str(sum(visual_response_stats_long == 1)) ' Visual Response 5 second window'])
disp([num2str(sum(visual_response_stats_short == 1 | visual_response_stats_long == 1)) ' Visual Response short or long'])
disp([num2str(sum(visual_response_stats_short == 0 & visual_response_stats_long == 1)) ' Visual Response long but not short'])
disp('-----------------------------------------------')
disp([num2str(sum((visual_response_stats_short == 1 | visual_response_stats_long == 1) & spatialness == 1)) ' Place cells are Visual Response short or long'])
disp([num2str(sum((visual_response_stats_short == 1) & spatialness == 1 & ~isnan(all_list_peak_times))) ' View Cells show Visual Response within 1 second'])
disp([num2str(sum((visual_response_stats_long == 1) & spatialness == 1 & ~ isnan(all_list_peak_times))) ' View Cells show Visual Response within 5 second'])
disp('-----------------------------------------------')
disp([num2str(sum(memory_long == 1 | memory_short == 1)) ' Show a Significant Memory Response'])
disp([num2str(sum((memory_long == 1 | memory_short == 1) & (spatialness == 1))) ' Place Cells Show a Significant Memory Response'])

%%
%---Copy Relevant Figures to Summary Directory---%
for unit = 1:length(all_unit_names)
    sub_dir1 = 'Visual Response\';
    
    %---Visual responsive neurons according to skaggs and correlation---%
    name1 = [all_unit_names{unit} '_Image_Visual Response.png'];
    if visual_response_stats_short(unit) == 1 || visual_response_stats_long(unit) == 1 %visually responsive
        if spatialness(unit) == 1 %place cell
            if ~exist([summary_directory 'Place\'],'dir')
                mkdir([summary_directory 'Place\']);
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Place\' name1])
        else %non place cell
            if ~exist([summary_directory 'Non Place\'],'dir')
                mkdir([summary_directory 'Non Place\']);
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Non Place\' name1])
        end
    end
    
    %---Visual responsive neurons only according to sliding window---%
    if (visual_response_stats_short(unit) == 0 && visual_response_stats_long(unit) == 0) ...
            && (visual_response_stats_short_sliding(unit) == 1 || visual_response_stats_long_sliding(unit) == 1)
        if ~exist([summary_directory 'Sliding Window Only\'],'dir')
            mkdir([summary_directory 'Sliding Window Only\']);
        end
        copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
            [summary_directory 'Sliding Window Only\' name1])
    end
    
    %---Fixation on Cross Response---%
    name2 = [all_unit_names{unit} '_Visual Response_CrossHair.png'];
    if fixation_on_cross_status(unit) == 1
        if ~exist([summary_directory 'Cross Response\'],'dir')
            mkdir([summary_directory 'Cross Response\']);
        end
        copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name2],...
            [summary_directory 'Cross Response\' name2])
    end
    
    %---Image Off Response---%
    if image_off_status(unit) == 1
        if ~exist([summary_directory 'Image Off Response\'],'dir')
            mkdir([summary_directory 'Image Off Response\']);
        end
        copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name2],...
            [summary_directory 'Image Off Response\' name2])
    end
    
    %---Memory Response---%
    name2 = [all_unit_names{unit} '_Image_Visual Response_Memory.png'];
    if memory_short(unit) == 1 %memory response within 1 second window
        if ~exist([summary_directory 'Memory\'],'dir')
            mkdir([summary_directory 'Memory\']);
        end
        copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name2],...
            [summary_directory 'Memory\' name2])
    end
    if memory_long(unit) == 1 %memory response within 1 second window
        if ~exist([summary_directory 'Memory\'],'dir')
            mkdir([summary_directory 'Memory\']);
        end
        copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name2],...
            [summary_directory 'Memory\' name2])
    end
    
end
%% Peri-Event Time Cell Plots for All Significant Cells Aligned to Max
figure
%---Time Cell Plot Fixation on Cross Hair---%
subplot(2,2,1)
[m,i] = max(cross_fixation_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin3-1],[1:size(cross_fixation_firing_rates,1)],cross_fixation_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(cross_fixation_firing_rates,1)],'w--');
hold off
xlabel('Time from Fixation  (ms)')
ylabel('Neuron #')
title('Crosshair Fixation')
xlim([-200 500])


%---Time Cell Plot Image On 1 second window---%
subplot(2,2,2)
[m,i] = max(image_on_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin2-1],[1:size(image_on_firing_rates,1)],image_on_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_on_firing_rates,1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Image Onset 1 s')
xlim([-200 1000])


%---Time Cell Plot Image On 5 second window---%
subplot(2,2,3)
[m,i] = max(long_image_on_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin1:twin4-1],[1:size(long_image_on_firing_rates,1)],long_image_on_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(long_image_on_firing_rates,1)],'w--');
hold off
xlabel('Time from Image On (ms)')
ylabel('Neuron #')
title('Image Onset 5 s')
xlim([-200 5000])

%---Time Cell Plot Image Off---%
subplot(2,2,4)
[m,i] = max(image_off_firing_rates,[],2);
[mm,ii] = sort(i);
imagesc([-twin3:twin3-1],[1:size(image_off_firing_rates,1)],image_off_firing_rates(ii,:))
colormap('jet')
hold on
plot([0 0],[1 size(image_off_firing_rates,1)],'w--');
hold off
xlabel('Time from Image Off (ms)')
ylabel('Neuron #')
title('Image Offset')
xlim([-500 500])