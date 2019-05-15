function Visual_Response_Memory(data_dir,figure_dir,session_data)
% modified from Visual_Response_AnalysisV2 to just do the memory component
% of the Visual response since both now require shuffling and may want to
% modify one without the other. Modified 12/27/16 SDK
%
% code rechecked for bugs January 15, 2017 by SDK

Fs = 1000;
task = 'ListSQ';
figure_dir = [figure_dir 'Visual Response\'];
numshuffs = 10000; %number of shuffles to do for bootstrapping
%want 10x more since want 2 tails @ p < 0.01 and at least 50 comparisons preferablly more

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task_file = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end

if ~exist([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat'],'file')
    return
end
load([data_dir task_file(1:8) '-ListSQ-Visual_Response_results.mat']);
num_units = size(unit_stats,2);
if num_units == 0
    return %if no units then return
end

%have to re-define here for parfor
%larger smoothing factor for longer window since want to remove "high"
%frequenies from analysis
smval = smval; %#ok
smval2 = smval2; %#ok


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Determine if Neuron Show Novel/Repeat Effect---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nov_rep_curve_short = cell(2,num_units);% for short 1 second window
sig_short = cell(1,num_units); %short sig times of p < 0.01
nov_rep_curve_long = cell(2,num_units); %for long 5 second window
sig_long = cell(1,num_units);  %long sig times of p < 0.01
for unit = 1:num_units
    %if neuron is showing stable response
    if ((epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)) || ...
            ((epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95))
        
        nov_rep_ind = find(~isnan(nvr{unit})); %can have NaNs which are meaningful here since unpaired novel/repeat images
        novreps = nvr{unit}; %all novel/repeat pair indeces
        novreps = novreps(nov_rep_ind); %remove unpaired trials
        
        %for short window, smooth first since saves computational times
        all_curves_short = NaN(numshuffs,twin1+twin2);
        [~,short_curves] = nandens(time_lock_firing{unit,3},smval,'gauss',Fs,'nanflt');%short window smoothed trial-by-trial
        short_curves = short_curves(~isnan(nvr{unit}),:); %remove unpaired trials
        
        %for long window, smooth first since saves computational times
        all_curves_long = NaN(numshuffs,twin1+twin4);
        [~,long_curves] = nandens(time_lock_firing{unit,5},smval2,'gauss',Fs,'nanflt');%short window smoothed trial-by-trial
        long_curves = long_curves(~isnan(nvr{unit}),:); %remove unpaired trials
        
        parfor shuff = 1:numshuffs
            ind = randperm(length(nov_rep_ind));
            
            %for short window
            nov_curve = nanmean(short_curves(novreps(ind) == 1,:)); %shuffled novel curve
            rep_curve = nanmean(short_curves(novreps(ind) == 2,:)); %shuffled repeat curve
            all_curves_short(shuff,:) = nov_curve-rep_curve; %shuffled difference
            
            %for long window
            nov_curve = nanmean(long_curves(novreps(ind) == 1,:)); %shuffled novel curve
            rep_curve = nanmean(long_curves(novreps(ind) == 2,:)); %shuffled repeat curve
            all_curves_long(shuff,:) = nov_curve-rep_curve; %shuffled difference
        end
        
        %---for short window---%
        nov_curve = nanmean(short_curves(novreps == 1,:)); %observed novel curve
        rep_curve = nanmean(short_curves(novreps == 2,:)); %observed repeat curve
        observed_diff = nov_curve-rep_curve; %observed_diff
        [~,sig_short{unit}] = cluster_level_statistic(observed_diff,all_curves_short,2,smval); %multiple comparision corrected significant indeces
        
        %---for long window---%
        nov_curve = nanmean(long_curves(novreps == 1,:)); %observed novel curve
        rep_curve = nanmean(long_curves(novreps == 2,:)); %observed repeat curve
        observed_diff = nov_curve-rep_curve; %observed_diff
        [~,sig_long{unit}] = cluster_level_statistic(observed_diff,all_curves_long,2,smval2); %multiple comparision corrected significant indeces
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot and Save Figures of Results---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t12 = -twin1:twin2-1;
t14 = -twin1:twin4-1;
unit_names = unit_stats(1,:);
for unit = 1:num_units
    if ((epoch_data.rate_prctile(unit,5) > 95) && (epoch_data.temporalstability_prctile(unit,5) > 95)) || ...
            ((epoch_data.rate_prctile(unit,3) > 95) && (epoch_data.temporalstability_prctile(unit,3) > 95))
        
        ylims = NaN(2,2);
        figure
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Short 1 second Window---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %---Plot raster for image on-short window---%
        subplot(2,3,1)
        [trial,time] = find(time_lock_firing{unit,3}(nvr{unit} == 1 | nvr{unit}== 2,:) == 1); %all pairs
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin2])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Plot Novel/Repeat raster for image on-short window---%
        subplot(2,3,2)
        hold on
        [trial,time] = find(time_lock_firing{unit,3}(nvr{unit} == 1,:) == 1); %novel
        if ~isempty(trial)
            plot(time-twin1,trial,'.b')
            xlim([-twin1 twin2])
            ylim([0 max(trial)+1]);
            b4 = max(trial);
        else
            b4 = 0;
        end
        [trial,time] = find(time_lock_firing{unit,3}(nvr{unit} == 2,:) == 1); %repeat
        if ~isempty(trial)
            trial = trial+b4;
            plot(time-twin1,trial,'.r')
            xlim([-twin1 twin2])
            ylim([0 max(trial)+1]);
        end
        hold off
        ylabel('Image #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Plot Firing rate locked to Image On-short window---%
        subplot(2,3,3)
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Long 5 Second Window---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Plot Raster for Long Window Image On---%
        subplot(2,3,4)
        [trial,time] = find(time_lock_firing{unit,5}(nvr{unit} == 1 | nvr{unit}== 2,:) == 1);%all pairs
        if ~isempty(trial)
            plot(time-twin1,trial,'.k')
            xlim([-twin1 twin4])
            ylim([0 max(trial)+1]);
        end
        ylabel('Trial #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Plot Novel/Repeat raster for image on-long window---%
        subplot(2,3,5)
        hold on
        [trial,time] = find(time_lock_firing{unit,5}(nvr{unit} == 1,:) == 1); %novel
        if ~isempty(trial)
            plot(time-twin1,trial,'.b')
            xlim([-twin1 twin4])
            ylim([0 max(trial)+1]);
            b4 = max(trial);
        else
            b4 = 0;
        end
        [trial,time] = find(time_lock_firing{unit,5}(nvr{unit} == 2,:) == 1); %repeat
        if ~isempty(trial)
            trial = trial+b4;
            plot(time-twin1,trial,'.r')
            xlim([-twin1 twin4])
            ylim([0 max(trial)+1]);
        end
        hold off
        ylabel('Image #')
        xlabel('Time from Image On (ms)')
        box off
        
        %---Plot Firing rate locked to Long Image On---%
        subplot(2,3,6)
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
        
        
        %---Set plots to same scale---%
        ymax = max(ylims(:,2));
        ymin = min(ylims(:,1));
        ymin(ymin < 0) = 0;
        
        subplot(2,3,3)
        ylim([ymin ymax])
        subplot(2,3,6)
        ylim([ymin ymax])
        
        %---Plot Significant Time Points---%
        subplot(2,3,3)
        hold on
        gaps = findgaps(find(sig_short{unit}));
        if ~isempty(gaps)
            for g = 1:size(gaps,1)
                gp = gaps(g,:);
                gp(gp == 0) = [];
                h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    [ymin ymin ymax ymax ymin],'k');
                uistack(h,'down')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
        end
        hold off
        
        subplot(2,3,6)
        hold on
        gaps = findgaps(find(sig_long{unit}));
        if ~isempty(gaps)
            for g = 1:size(gaps,1)
                gp = gaps(g,:);
                gp(gp == 0) = [];
                h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
                    [ymin ymin ymax ymax ymin],'k');
                uistack(h,'down')
                set(h,'facealpha',.25,'EdgeColor','None')
            end
        end
        hold off
        
        subtitle(['Visual Response/Memory ' task_file(1:8) '_' unit_names{unit} ', n_ = ' num2str(sum(nvr{unit} == 1))]);
        save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names{unit} '_Image_Visual Response_Memory']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Finally save all the data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-ListSQ-Visual_Response_Memory_results.mat'],...
    'smval','smval2','sig_short','sig_long','numshuffs','unit_names')
disp(['Memory Analyis for ' task_file(1:8) ' saved']);
end