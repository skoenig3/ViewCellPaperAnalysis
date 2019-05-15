function valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task)
%written by Seth Konig. Since use so much put in own file
%for each unit determine the which trials are valid...if the unit is not
%well isolated for min_blks then throw out unit
%
% Code rechecked ListSQ section for bugs October 17, 2016 SDK

if strcmpi(task,'cvtnew')
    trials_per_block = 100; %approximately, varies slightly depending on path length
    
    num_trials = length(cfg.trl);
    %NaNs are for start and end trials otherwise cut
    valid_trials(1,isnan(valid_trials(1,:))) = 1;
    valid_trials(2,isnan(valid_trials(2,:))) = num_trials;
    
    for unit = 1:num_units
        start_end = valid_trials(:,unit);
        if isnan(start_end(1)) %NaNs meant first trial was usable
            start_end(1) = 1;
        end
        if isnan(start_end(2))%NaNs meant last trial was usable
            start_end(2) = length(cfg.trl);
        end
        start_end(start_end == 0) = 1; %set 0 to 1 for indexing during next line
        num_blks = floor((start_end(2)-start_end(1)+1)/trials_per_block);
        
        if num_blks < min_blks
            valid_trials(:,unit) = NaN;
        end
    end
    
    
elseif strcmpi(task,'ListSQ')
    num_trials = length(cfg.trl);
    %NaNs are for start and end trials otherwise cut
    valid_trials(1,isnan(valid_trials(1,:))) = 1;
    valid_trials(2,isnan(valid_trials(2,:))) = num_trials;
    
    if str2double(task_file(3:8)) < 140805 %7 second viewing with fewere sequence trials
       trials_per_block = 80;%for rest of session: includes ~16 novel/repeat images + sequence trials
    else
       trials_per_block =  96;%for rest of session: includes ~16 novel/repeat images + sequence trials
    end
    
    for unit = 1:num_units
        start_end = valid_trials(:,unit);
        if isnan(start_end(1)) %NaNs meant first trial was usable
            start_end(1) = 1;
        end
        if isnan(start_end(2))%NaNs meant last trial was usable
            start_end(2) = length(cfg.trl);
        end
        start_end(start_end == 0) = 1; %set 0 to 1 for indexing during next line
        min_trial = cfg.trl(start_end(1)).cnd-1000; %get the first condition number
        max_trial = cfg.trl(start_end(2)).cnd-1000; %get the last condition number
        
        if min_trial < 22 %includes familizarization block for sequences
            min_trial = 22; %remove then count from there because doesn't include pictures
        end
        num_blks = floor((max_trial-min_trial+1)/trials_per_block);
        
        if num_blks < min_blks
            valid_trials(:,unit) = NaN;
        end
    end
else
    error('Could not identify task')
end
