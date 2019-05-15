function spatial_analysis_plotsV2(figure_dir,task_file,position,spike_times,spatial_info,...
    task_type,unit_names,smval,imageX,imageY,trial_type,Fs,peak_firing_rate)
% written by Seth Konig August, 2014
% generates and save plots for spatial spike analysis
% % updated SDK 1/11/17 to handlde new format and partial session data for
% vaild trials only. Updated CVTNEW section on 1/19/16
%
% Inputs:
%   1) figure_dir: directory where to save figures
%   2) task_file: name of the preprocessed_file
%   3) position: of attention by trial x in odd rows y in even rows
%   4) spike_times: aligned to eye by col and trial in row
%   6) spatial_info: spatial information
%   5) task_type: type of plot to create by task and/or subtask...
%       a) 'List_spatial': plot firing rate over space
%       b) 'List_spatial_time_shifted': time shift variations of 'List_spatial'
%       c) 'cvtnew_spatial': plot firing rate over space
%   6) unit_names...
%       a) unit_names.name: name of unit e.g. sig001a
%       b) unit_names.multiunit: 1 for multiunit 0 for single unit
%   7) smval: smoothing parameters of varying lengths depending on task
%   8) imageX: horizontal size of the image
%   9) imageY: vertical size of the image
%   10) trial_type: for List task only e.g. novel vs repeat
%   11) Sampling rate
%
% Outputs:
%   1) saved figures into figure_dir
%
% Code rechecked ListSQ section for bugs October 18-19, 2016 SDK

figure_dir = [figure_dir 'Spatial Analysis\'];
num_units = length(position);

switch task_type
    case 'List_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        H = define_spatial_filter(filter_width);
        
        for unit = 1:num_units
            if ~isnan(peak_firing_rate(1,unit))
                clims = NaN(2,3);
                num_trials_str = [' n_ = ' num2str(size(spike_times{unit},1))];
                
                figure
                
                %---draw raw spikes on all images all locations---%
                make_spike_jittered_plot(position{unit},spike_times{unit},[2 3],1)
                title(['Raw: All images,' num_trials_str])
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                %---plot spike locations color coded by 1st half vs second half of session/by time---%
                make_spike_jittered_colored_plot(position{unit},spike_times{unit},[2 3],2)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                title(sprintf(['Halves: \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(unit),2) ...
                    ' ' num2str(spatial_info.spatialstability_halves_prctile(unit),3) '%%']))
                
                %---plot spike locations color code by even and odd trials---%
                make_spike_jittered_colored_plot_even_odd(position{unit},spike_times{unit},[2 3],3)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                
                title(sprintf(['Even/Odd: \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(unit),2) ...
                    ' ' num2str(spatial_info.spatialstability_even_odd_prctile(unit),3) '%%']))
                
                
                %---Plot Rate Map for All Images---%
                ratemap = get_firing_rate_map({position{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');
                
                subplot(2,3,4)
                h = imagesc(ratemap);
                %                 hold on
                %                 gry = 0.5*ones(size(ratemap,1),size(ratemap,2),3);
                %                 h = imshow(gry);
                %                 hold off
                %                 set(h,'alphadata',isnan(ratemap));
                set(h,'alphadata',~isnan(ratemap));
                axis off
                axis equal
                
                clims(:,1) = caxis;
                max_fr = prctile(ratemap(:),97.5); % the ~97.5%-tile
                clims(2,1) = max_fr;
                
                title(sprintf(['All images, peak rate = ' num2str(max_fr,3) ...
                    ' Hz \n Bit/s = ' num2str(spatial_info.rate(unit),3) ' ' ...
                    num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%']))
                
                %---plot firing rate map for novel images---%
                nov_ratemap = get_firing_rate_map({select_eyepos(position{unit},trial_type{unit} ==1),...
                    spike_times{unit}(trial_type{unit} ==1,:)},imageX,imageY,binsize,H,Fs,'all');
                
                subplot(2,3,5)
                h = imagesc(nov_ratemap);
                set(h,'alphadata',~isnan(nov_ratemap));
                axis off
                axis equal
                
                
                clims(:,2) = caxis;
                max_fr = prctile(nov_ratemap(:),97.5); % the ~97.5%-tile
                clims(2,2) = max_fr;
                
                title(sprintf(['Novel images, peak rate = '  num2str(max_fr,3) ' Hz']))
                
                
                %---plot firing rate map for repeat images---%
                rep_ratemap = get_firing_rate_map({select_eyepos(position{unit},trial_type{unit} ==2),...
                    spike_times{unit}(trial_type{unit} ==2,:)},imageX,imageY,binsize,H,Fs,'all');
                
                subplot(2,3,6)
                h = imagesc(rep_ratemap);
                set(h,'alphadata',~isnan(rep_ratemap));
                axis off
                axis equal
                
                
                clims(:,3) = caxis;
                max_fr = prctile(rep_ratemap(:),97.5); % the ~97.5%-tile
                clims(2,3) = max_fr;
                
                title(sprintf(['Repeat images, peak rate = '  num2str(max_fr,3) ' Hz']))
                
                minc = nanmin(clims(1,:));
                maxc = nanmax(clims(2,:));
                diffc = maxc-minc;
                minc = minc - 0.15*diffc;
                for sp = 4:6
                    subplot(2,3,sp)
                    %colormap(viridis)
                    colormap('jet')
                    caxis([minc maxc])
                    
                    if sp == 4;
                        colorbar
                    end
                    
                    c = colormap;
                    c(1,:) = [1 1 1];%turn nans in to white pixels
                    colormap(c);
                end
                
                %label as putative single or multi-unit
                if unit_names.multiunit(unit)
                    multi_str = 'Multiunit ';
                else
                    multi_str = ' ';
                end
                
                %label as putative excitatory/inhibitory neuron based on
                %whole session average firing rate
                if unit_names.putative_EI(unit) == 1 %putative excitatory neuron
                    avg_firing_rate_str = [', avg FR = ' num2str(unit_names.avg_firing_rate(unit),2) ...
                        ' Hz, putative Excitatory Neuron'];
                elseif unit_names.putative_EI(unit) == 2 %putative inhibitory neuron
                    avg_firing_rate_str = [', avg FR = ' num2str(unit_names.avg_firing_rate(unit),2) ...
                        ' Hz, putative Inhibitory Neuron'];
                end
                
                subtitle([multi_str task_file(1:8) ' ' multi_str unit_names.name{unit} avg_firing_rate_str]);
                
                if unit_names.multiunit(unit)
                    save_and_close_fig([figure_dir '\MultiUnit\'],[task_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
                else
                    save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_List_spatial_analysis']);
                end
            end
        end
        
    case 'cvtnew_spatial'
        binsize = smval(1); %number of pixels per spatial bing
        filter_width = smval(2); %std of 2D gaussian smoothing filter
        H = define_spatial_filter(filter_width);
        H2 = fspecial('gaussian',6*1+1,1);%for filtering eye data on crosshair
        
        
        for unit = 1:num_units
            
            if ~isempty(position{unit})
                
                correct_trials = unit_names.correct_trials{unit};
                
                %remove data for incorrect trials since may not have been paying attention
                dtpsx =position{unit}(1:2:end,:);
                dtpsx = dtpsx(correct_trials == 1,:);
                dtpsy = position{unit}(2:2:end,:);
                dtpsy = dtpsy(correct_trials == 1,:);
                
                dtps = NaN(2*size(dtpsx,1),size(dtpsx,2));
                dtps(1:2:end,:) = dtpsx;
                dtps(2:2:end,:) = dtpsy;
                
                position{unit} = dtps;
                spike_times{unit} = spike_times{unit}(correct_trials == 1,:);
                
                
                figure
                
                %spike postion overlaying dot position
                make_spike_jittered_plot(position{unit},spike_times{unit},[3 3],1)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                xlim([100 700])
                ylim([0 600])
                axis off
                axis square
                
                %firing ratemap
                [ratemap,~] = get_firing_rate_map_cvtnew({position{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');
                maxfr = prctile(ratemap(:),97.5);
                subplot(3,3,4)
                h = imagesc(ratemap);
                set(h,'alphadata',~isnan(ratemap));
                axis off
                axis square
                colorbar
                colormap('jet')
                clim = caxis;
                caxis([clim(1) maxfr])
                
                title_str = ['peak rate = ' num2str(maxfr,3) ' Hz'];
                if spatial_info.shuffled_rate_prctile(unit) > 90;
                    title_str = [title_str ' \n Bits ' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%'];
                end
                if spatial_info.spatialstability_halves_prctile(unit) > 90;
                    title_str = [title_str  '\n r = ' num2str(spatial_info.spatialstability_halves(unit),2) ' ' ...
                        num2str(spatial_info.spatialstability_halves_prctile(unit),2) '%%'];
                end
                title(sprintf(title_str))
                
                path_number = unit_names.path_number{unit}(correct_trials == 1);
                %spike position overlaying dot position for 1st path recorded
                this_path = find(path_number == 1);
                %num_trials = floor(size(spike_times{unit},1)/2);
                make_spike_jittered_plot(position{unit}((2*this_path(1)-1):this_path(end)*2,:),spike_times{unit}(this_path,:),[3 3],2)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                xlim([100 700])
                ylim([0 600])
                axis off
                axis square
                title('1st Path')
                
                %spike position overlaying dot position for 2nd path recorded
                this_path = find(path_number == 2);
                make_spike_jittered_plot(position{unit}((2*this_path(1)-1):this_path(end)*2,:),spike_times{unit}(this_path,:),[3 3],5)
                set(gca,'Xcolor','w')
                set(gca,'Ycolor','w')
                xlim([100 700])
                ylim([0 600])
                axis off
                axis square
                title('2nd path')
                
                num_not_nans = sum(~isnan(spike_times{unit}));
                indeces = find(num_not_nans > 25);
                
                %firing rate curve aligned to do on
                subplot(3,3,3)
                dofill2(1:length(indeces),spike_times{unit}(:,indeces),'black',1,60);
                xlim([0 2000])
                yl = ylim;
                if yl(1) < 0
                    ylim([0 yl(2)])
                end
                xlabel('Time from Dot On (ms)')
                ylabel('Firing Rate (Hz)')
                
                %raster aligned to dot on
                subplot(3,3,6)
                [trial,time] = find(spike_times{unit}(:,indeces) == 1);
                if ~isempty(trial)
                    plot(time,trial,'.k')
                    xlim([0 2000])
                    ylim([0 max(trial)+1]);
                end
                ylabel('Trial #')
                xlabel('Time from Dot On (ms)')
                box off
                
                subplot(3,3,7)
                imagesc(log10(imfilter(unit_names.eyepos{1,unit},H2)))
                hold on
                rectangle('Position',[400-54/2 300-54/2 54 54],'EdgeColor',[1 1 1],'LineWidth',2)
                hold off
                xlim([350 450])
                ylim([250 350])
                axis off
                axis square
                title('Log_{10} All Eye Position')
                
                subplot(3,3,8)
                imagesc(log10(imfilter(unit_names.eyepos{2,unit},H2)))
                hold on
                rectangle('Position',[400-54/2 300-54/2 54 54],'EdgeColor',[1 1 1],'LineWidth',2)
                hold off
                xlim([350 450])
                ylim([250 350])
                axis off
                axis square
                title('Log_{10} Eye Pos. When Spike')
                
                subplot(3,3,9)
                plot(-1000:1000,spatial_info.observed_STA{unit},'k')
                hold on
                plot(-1000:1000,spatial_info.shuffled_95_direction_STA{unit},'r')
                hold off
                xlabel('Time')
                ylabel('MRL')
                title('Direction STA')
                box off
                
                
                if unit_names.multiunit(unit)
                    subtitle(['Multiunit ' task_file(1:end-11) ' ' unit_names.name{unit} ]);
                else
                    subtitle([task_file(1:end-11) ' ' unit_names.name{unit}]);
                end
                save_and_close_fig(figure_dir,[task_file(1:end-11) '_' unit_names.name{unit} '_cvtnew_spatial_analysis']);
            end
        end
end
end