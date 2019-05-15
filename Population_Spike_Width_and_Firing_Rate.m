% Code to calculate spike width and firing to determine interneuron
%vs pyramidal neuron. Single units only
clar
task = 'ListSQ';

mean_firing_rate = [];
peak_valley_width = [];
half_max_width = [];
median_ISI = [];
percent_ISI_25 = [];
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
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
            if multiunit(unit) == 1 || isnan(multiunit(unit)) || isnan(valid_trials(1,unit))
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
            
            avg_waveform = mean(waveforms{1,unit},2);
            [mx,mxi] = max(avg_waveform);
            [mn,mni] = min(avg_waveform);
            if max(avg_waveform) > abs(2*min(avg_waveform)) %call it  inverted
                half_max = 1/2*mx;
                above_half = find(avg_waveform > half_max);
                if avg_waveform(1) > half_max  %started off large
                    disp('now');
                else
                     disp('now');
                end
                
            else %"normal" unit
                half_min = 1/2*mn; 
                below_half = find(avg_waveform < half_min);
                half_max_width = [half_max_width [below_half(end)-below_half(1)+1]];                
            end     
            peak_valley_width = [peak_valley_width abs(mxi-mni)];
       
            
            ISIS(ISIS > 500) = [];
            
            
          
            
            median_ISI = [median_ISI median(ISIS)];
            percent_ISI_25 = [percent_ISI_25 sum(ISIS < 25)/length(ISIS)];
            
            time_count = time_count/1000;
            mean_firing_rate = [mean_firing_rate spike_count/time_count];
            
            %             if spike_count/time_count > 10 && sum(ISIS < 25)/length(ISIS) > 0.6
            %                 disp('now')
            %             end
        end
        
        
    end
end
% 25 micro sec per sample
%%
figure
subplot(1,2,1)
plot(mean_firing_rate,25*half_max_width,'k.')
xlabel('Mean Firing Rate (Hz)')
ylabel('FWHM (us)')

subplot(1,2,2)
plot(mean_firing_rate,25*peak_valley_width,'k.')
xlabel('Mean Firing Rate (Hz)')
ylabel('Peak-2-Valley Width (us)')

%%
figure
mfr = mean_firing_rate;
mfr(mfr > 30) = 30;
hist(mfr,50)
xlabel('Mean Firing Rate (Hz)')
ylabel('Count')
%%
IDK = kmeans([mfr' half_max_width' peak_valley_width'],2);
IDX = kmeans([mfr'],2);
IDXX = kmeans([mfr' peak_valley_width'],2);
%%
IDKK =  kmeans([mean_firing_rate' half_max_width' percent_ISI_25'],3);

%%
clr = 'rgb';
figure
hold on
subplot(1,2,1)
hold on
for c = 1:3
    plot(mean_firing_rate(IDKK == c),25*half_max_width(IDKK == c),['.' clr(c)])
end
xlabel('Mean Firing Rate (Hz)')
ylabel('FWHM (us)')


subplot(1,2,2)
hold on
for c = 1:3
    plot(mean_firing_rate(IDKK == c),percent_ISI_25(IDKK == c),['.' clr(c)])
end
xlabel('Mean Firing Rate (Hz)')
ylabel('% ISI < 25 ms')
%%
fast_units = mfr(IDX == 2);
PI2 = percent_ISI_25(IDX == 2);
mI = median_ISI(IDX == 2);
FWHM =half_max_width(IDX == 2);

plot3(fast_units,PI2,FWHM,'k.')
xlabel('Firing Rate (Hz)')
ylabel('% ISI < 25 ms')
zlabel('FWHM')
%%
figure
subplot(1,2,1)
plot(mean_firing_rate,25*half_max_width,'k.')
hold on
plot([0 8.5],[mean(25*half_max_width(IDX == 1)) mean(25*half_max_width(IDX == 1))],'r')
plot([8.5 45],[mean(25*half_max_width(IDX == 2)) mean(25*half_max_width(IDX == 2))],'b')
hold off
xlabel('Mean Firing Rate (Hz)')
ylabel('FWHM (us)')

subplot(1,2,2)
plot(mean_firing_rate,25*peak_valley_width,'k.')
hold on
plot([0 8.5],[mean(25*peak_valley_width(IDX == 1)) mean(25*peak_valley_width(IDX == 1))],'r')
plot([8.5 45],[mean(25*peak_valley_width(IDX == 2)) mean(25*peak_valley_width(IDX == 2))],'b')
hold off
xlabel('Mean Firing Rate (Hz)')
ylabel('Peak-2-Valley Width (us)')

%%
figure
plot(mean_firing_rate,percent_ISI_25,'k.')
hold on
plot([0 8.5],[mean(percent_ISI_25(IDX == 1)) mean(percent_ISI_25(IDX == 1))],'b')
plot([8.5 45],[mean(percent_ISI_25(IDX == 2)) mean(percent_ISI_25(IDX == 2))],'r')
plot(15.4,0.7965,'go','markersize',10)
xlabel('Firing Rate (Hz)')
ylabel('% ISI < 25 ms')

%%
figure
plot(mean_firing_rate,median_ISI,'k.')
xlabel('Firing Rate (Hz)')
ylabel('Median ISI (ms)')
