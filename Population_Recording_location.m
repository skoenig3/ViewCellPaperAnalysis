%get recording locations for view cells
clar

task = 'ListSQ';


all_unit_AP_location = NaN(1,347);
all_unit_subregion = cell(1,347);
unit_count = 1;

monkeys = {'Vivian','Tobii'};
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
    
    for session = 1:length(session_data)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{session},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        clear valid_trials
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'waveforms','cfg','valid_trials',...
            'data','item_file','cnd_file','hdr');
        %get important task specific information
        %grab unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        
        if ~exist('valid_trials') && num_units == 0
            continue
        elseif  ~exist('valid_trials') && num_units ~= 0
            error('Should have units')
        end
        
        num_trials = length(cfg.trl);
        %NaNs are for start and end trials otherwise cut
        valid_trials(1,isnan(valid_trials(1,:))) = 1;
        valid_trials(2,isnan(valid_trials(2,:))) = num_trials;
        
        if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
            minimum_trials_1 = 101; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
            minimum_trials_2 =  80;%for rest of session: includes 16 novel/repeat images + sequence trials
        else
            minimum_trials_1 = 117; %for session start minimum number of trials: includes fam block and 16 novel/repeat images + sequence trials
            minimum_trials_2 =  96;%for rest of session: includes 16 novel/repeat images + sequence trials
        end
        
        min_blks = 2;
        for unit = 1:num_units
            start_end = valid_trials(:,unit);
            if isnan(start_end(1))
                start_end(1) = 1;
            end
            if isnan(start_end(2))
                start_end(2) = length(cfg.trl);
            end
            start_end(start_end == 0) = 1;
            min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
            max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
            
            if min_trial < 22 %includes fam block
                min_trial = 22; %remove then count from there
            end
            num_blks = floor((max_trial-min_trial+1)/minimum_trials_2);
            
            if num_blks < min_blks
                valid_trials(:,unit) = NaN;
            end
        end
        
        for unit = 1:num_units
            if isnan(multiunit(unit)) || isnan(valid_trials(1,unit))
                continue
            end
            
            all_unit_AP_location(unit_count) = chamber_zero(1)+ session_data{session}.location(1);
            all_unit_subregion{unit_count} = session_data{session}.subregion;
            unit_count = unit_count + 1;
        end
    end
end
%%

figure
hist(all_unit_AP_location,15)
xlabel('Anterior-Posterior Location (mm)')
yalbel('Unit count')

%%
DG_count = sum(contains(all_unit_subregion,'DG') | contains(all_unit_subregion,'CA4'));
CA3_count = sum(contains(all_unit_subregion,'CA3') | contains(all_unit_subregion,'CA2'));
CA1_count = sum(contains(all_unit_subregion,'CA1') | contains(all_unit_subregion,'sub'));