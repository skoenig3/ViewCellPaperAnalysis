% Code below creates population summary for Significnat Saccade Amplitude Cells
% Written by Seth Konig modified from Population_Saccade_amplitude. on 1/31/2017
% Code does the following

clar %clear,clc

%where to store spatial analysis figure copies
summary_directory = 'D:\MATLAB\ViewCellPaperAnalysis\PopulationFigures\Saccade Amplitude\';
if ~isdir(summary_directory)
    mkdir(summary_directory)
end

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
image_on_twin = 500;%ms

%---Values All Fixations out2out and in2in---%
all_corrs = []; %all correlation between firing rate and saccade amplitude
all_corrs_pctiles = []; %observed shuffled percentile

%---Other values---%
spatialness = []; %1 for place cell, 0 for non place cell
all_unit_names = {};
all_monkeys = []; %1s and 2s for monkeys
all_sacamps = [];

monkeys = {'Vivian','Tobii'};
figure_dir = {};
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = 'P:\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Recording Files\';
        figure_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = 'P:\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
        figure_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess =1:length(session_data)
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
            load([data_dir  task_file(1:8) '-Eyemovement_Locked_List_results.mat'],'fixation_information');
        else
            if num_units ~= 0
                error('Where is this file')
            end
            continue
        end
        
        
        num_units = size(unit_stats,2);
        for unit = 1:num_units
            if ~isnan(amplitude_correlations(unit)) %if unit was processed
                
                %---Get Distribution of Saccade Amplitudes---%
                sac_amps = fixation_information{unit}(:,6); %fixation aligned firing rate
                sac_amps = sac_amps/24; %convert to dva
                
                %remove fixation within the first 500 ms of images on
                %amplitude is affected by time within image period
                time_from_image_on = fixation_information{unit}(:,4); %fixation aligned firing rate
                too_early = find(time_from_image_on < image_on_twin);
                sac_amps(too_early,:) = [];
                
                all_sacamps = [all_sacamps; sac_amps]; %store across sessions
                
                %---Values All Fixations out2out and in2in---%
                all_corrs = [all_corrs amplitude_correlations(unit)]; %all observed correlation
                all_corrs_pctiles = [all_corrs_pctiles amplitude_correlations_percentile(unit)];%observed shuffled percentile

                
                %---Other values---%
                all_unit_names = [all_unit_names {[task_file(1:end-11) '_' unit_stats{1,unit}]}];
                all_monkeys = [all_monkeys monkey];
                if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
                    spatialness = [spatialness 1]; %place cell
                else
                    spatialness = [spatialness 0]; %non place cell
                end
            end
        end
    end
end
%%
clc
disp([num2str(nansum(spatialness)) ' place cells'])
disp([num2str(sum(all_corrs_pctiles > 97.5)) ' amplitude modulated cells'])
disp('--------------------------------------------------------------')
disp([num2str(sum(all_corrs_pctiles > 97.5 & spatialness == 1)) ' amplitude modulated place cells'])
disp([num2str(sum(all_corrs_pctiles > 97.5 & spatialness == 0)) ' amplitude modulated non-place cells '])
disp('--------------------------------------------------------------')
disp([num2str(sum(~isnan(all_corrs_pctiles))) ' in total analyzed'])
disp([num2str(sum((~isnan(all_corrs_pctiles) & (spatialness == 1)))) ' place cells analyzed'])
disp([num2str(sum((~isnan(all_corrs_pctiles) & (spatialness == 0)))) ' non-place cells analyzed'])

%%
%---Copy Relevant Figures to Summary Directory---%
for unit = 1:length(all_unit_names)
    if all_corrs_pctiles(unit) > 97.5
        sub_dir1 = 'Saccade Direction and Amplitude\';
        name1 = [all_unit_names{unit} '_Saccade_Amplitude_Analysis.png'];
        if spatialness(unit) == 1 %place cell
            if ~exist([summary_directory 'Place\'],'dir')
               mkdir([summary_directory 'Place\']); 
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Place\' name1])
        elseif spatialness(unit) == 0 %non place cell
            if ~exist([summary_directory 'Non Place\'],'dir')
               mkdir([summary_directory 'Non Place\']); 
            end
            copyfile([figure_dir{all_monkeys(unit)} sub_dir1 name1],...
                [summary_directory 'Non Place\' name1])
        end
    end
end
