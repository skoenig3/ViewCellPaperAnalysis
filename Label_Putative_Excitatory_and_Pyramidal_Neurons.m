%code written by Seth Konig August 22, 2016
%code labels neurons as putative interneurons or exctiatory neurons based
%on mean firing rate across the whole time that the neuron is stable. MU
%and SU are labeled the same though MU will be less accurate. However,
%internurons in MU will make the MU look inhibitory anway while multiple
%excitatory neurons may not. We could change MU criterion based on number
%of clusters in PCA or something else.


clar
task = 'ListSQ';

firing_rate_threshold = 8.5; %based on 2 clusters on SU with 2+ blocks of stable activity
%in both Vivian and Tobii (similar for them individually). On average units
%with wider waveforms were slow but not as clean cut so only using firing
%rate

ISI_threshold = 0.75;% if % of ISI < 25 ms is > 75% and unit is above the
%firing rate_threshold then will still call excitatory neuron may just be
%bursty and very responsive.

min_blks = 1; %will will label neurons that have at least 1+ novel/repeat blocks in case we use them later

excitatory_inhibitory_counts = zeros(1,3);
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for session = 1:length(session_data);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
            = get_task_data(session_data{session},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        
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
        
        excitatory_inhibitory = NaN(1,num_units);%1 for excitatory 2 for inhibitory, 3 possible fast firing inhibitory neuron
        whole_session_mean_firing_rate = NaN(1,num_units);
        for unit = 1:num_units
            if isnan(valid_trials(1,unit))
               continue 
            end
            
            spike_count = 0;
            time_count = 0;
            ISIS = [];
            for t = 1:num_trials
                if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
                    if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                        spike_count = spike_count + nansum(data(unit).values{t});
                        time_count = time_count +length(data(unit).values{t});
                        spk = find(data(unit).values{t});
                        ISIS = [ISIS diff(spk)];
                    elseif any(cfg.trl(t).allval == 23)  && any(cfg.trl(t).allval == 3) %rewarded sequence trials
                        spike_count = spike_count + nansum(data(unit).values{t});
                        time_count = time_count +length(data(unit).values{t});
                        spk = find(data(unit).values{t});
                        ISIS = [ISIS diff(spk)];
                    end
                end
            end
            
            time_count = time_count/1000;
            mean_firing_rate =  spike_count/time_count;
            whole_session_mean_firing_rate(unit) = mean_firing_rate;
            
            ISIS(ISIS > 500) = [];
            percent_ISI_25 = sum(ISIS < 25)/length(ISIS);
            
            if mean_firing_rate < firing_rate_threshold %then putative excitatory neuron
                 excitatory_inhibitory(unit) = 1;
                 excitatory_inhibitory_counts(1) = excitatory_inhibitory_counts(1)+1;
            elseif  percent_ISI_25 > ISI_threshold %possibly fast firing excitatory
                excitatory_inhibitory(unit) = 3;
                excitatory_inhibitory_counts(3) = excitatory_inhibitory_counts(3)+1;
            else %putative inhibitory neuron
                excitatory_inhibitory(unit) = 2;
                excitatory_inhibitory_counts(2) = excitatory_inhibitory_counts(2)+1;
            end
        end
        save([data_dir task_file(1:end-11) '-preprocessed.mat'],'-append',...
            'excitatory_inhibitory','whole_session_mean_firing_rate'); %append to preprocesed file
    end
end
