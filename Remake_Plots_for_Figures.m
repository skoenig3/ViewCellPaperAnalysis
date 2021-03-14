% %% recreate plots so can pull them into illustrator
%
%% Just List


% %plots fixation aligned rasters so that can export to eps
clar
data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';

task_file = 'TO160325';
unit_name = 'sig002a';

load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

fix_locked_firing = list_fixation_locked_firing{this_unit};
fix_in_out = in_out{this_unit};
t = -twin1:twin2-1;

num_in = sum(fix_in_out == 1);
num_out = sum(fix_in_out == 4);
downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number
figure

%---Fixations in->out vs out->out---%
out_matrix = fix_locked_firing(fix_in_out == 4,:);
out_matrix = out_matrix(1:downsample:end,:);
[trial,time] = find(out_matrix == 1);
plot(time-twin1,(trial),'.b')
hold on
if ~isempty(trial)
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1);
trial = trial+b4;
plot(time-twin1,(trial),'.r')
if ~isempty(trial)
    ylim([0 max(trial)])
else
    ylim([0 b4])
end
box off
plot([0 0],[0 max(trial)+1],'k--')
ylabel('Occurence')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
xlabel('Time from Fixation Start (ms)')
title('Fixation Aligned Rasters')
axis square

%---Firing Rate Curves---%
figure
hold on
[~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
%     plot(t,list_95_curve{1,unit},'k','linewidth',2);%95% confidence interval
[pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
locs(pks < 0.66*max(y_list)) = [];
pks(pks < 0.66*max(y_list)) = [];
plot(locs-twin1,pks,'*k')
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
end
plot([0 0],[yl(1) yl(2)],'k--')
% gaps = findgaps(sig_ind);
% if ~isempty(gaps)
%     for g = 1:size(gaps,1)
%         gp = gaps(g,:);
%         gp(gp == 0) = [];
%         if length(gp) > 40
%             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
%                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
%             uistack(h,'down')
%             set(h,'facealpha',.25,'EdgeColor','None')
%         end
%     end
% end
xlim([-twin1 twin2]);
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
legend('out->in','out->out','Location','NorthWest')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
% title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
%     ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
axis square


%---Normalized Firing Rate Curve---%
%for plot2b
clear var area %exist in store data 

%normalize out2in firing rate curve
in_curve = y_list;
in_curve = in_curve-mean(in_curve(1:twin1));
in_curve = in_curve/max(in_curve);

figure
area(-twin1:twin2-1,in_curve)
hold on
plot([-twin1 twin2],[0 0],'k')
yls = ylim;
plot([0 0],[yls(1) yls(2)],'k--')
hold off

%%
%% List vs Sequence Analysis
%
clar
data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
%data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Resorted\';

task_file = 'TO160318';
unit_name = 'sig003a';

% load([data_dir  task_file '-preprocessed.mat'],'cfg','item_file','cnd_file');
load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
load([data_dir  task_file(1:8) '_3-spatial_analysis_results.mat']);
load([data_dir  task_file(1:8) '_3-preprocessed.mat'],'item_file','cnd_file');
[~,~,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);

H = define_spatial_filter(filter_width);


this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end
imageX = 800;
imageY = 600;

firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
maxfr = prctile(firing_rate_map(:),97.5);

figure
%---plot firing rate map for all images---%
h = imagesc(firing_rate_map);
set(h,'alphadata',~isnan(firing_rate_map));
title('All images')
axis off
axis equal
colorbar
colormap('jet')
clim = caxis;
caxis([clim(1) maxfr])

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

fix_locked_firing = list_fixation_locked_firing{this_unit};
fix_in_out = in_out{this_unit};
t = -twin1:twin2-1;

num_in = sum(fix_in_out == 1);
num_out = sum(fix_in_out == 4);
downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number

place_field_matrix = all_place_field_matrix{unit};
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

figure
hold on
[~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
%     plot(t,list_95_curve{1,unit},'k','linewidth',2);%95% confidence interval
[pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
locs(pks < 0.66*max(y_list)) = [];
pks(pks < 0.66*max(y_list)) = [];
plot(locs-twin1,pks,'*k')
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
end
plot([0 0],[yl(1) yl(2)],'k--')
%
xlim([-twin1 twin2]);
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
legend('out->in','out->out','Location','NorthWest')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
% title(sprintf(['n_{out->in} = ' num2str(sum(fix_in_out == 1))...
%     ', n_{out->out} = ' num2str(sum(fix_in_out == 4))]));
title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
axis square

%---Plot Place Field---%
colors = 'rg';
shapes = 'xo';

figure
imagesc(all_place_field_matrix{unit});
colormap('gray')
hold on
for c = 1:4
    for seq = 1:2
        if ~isnan(sequence_inside(seq,c))
            plot(sequence_locations{seq}(1,c),imageY-sequence_locations{seq}(2,c),[colors(seq) shapes(seq)],'markersize',16)
        end
    end
end
hold off
xlim([0 800])
ylim([0 240])
axis equal
axis off
title(sprintf(['Place Field Location \n Area = ' num2str(area(unit),2) '%%']));


%---Plot Firing Rate Curves for Suquence Trials---%
which_sequence = all_which_sequence{unit};
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

figure
hold on
dofill(t,fixation_firing(in_out_sequence{unit} == 1,:),'red',1,smval);
dofill(t,fixation_firing(in_out_sequence{unit} == 0,:),'blue',1,smval);
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
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



%---Plot Firing Rate Curves for Suquence vs Image Trials---%
figure
hold on
dofill(t(1:twin1+twin2),list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'black',1,smval);%list out-> in
dofill(t,fixation_firing(in_out_sequence{unit} == 1,:),'green',1,smval);%sequence
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

if stats_across_tasks(2,unit) > stats_across_tasks(4,unit)
    title(['Contextual Gain: ' num2str(100*(stats_across_tasks(2,unit)/stats_across_tasks(4,unit)),3) '%'])
else
    title(['Contextual Gain: ' num2str(-100*(stats_across_tasks(4,unit)/stats_across_tasks(2,unit)),3) '%'])
end

%% Saccade Direction Stuff

%plots fixation aligned rasters so that can export to eps

clar
data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
%data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Resorted\';

task_file = 'TO151208_3';
unit_name = 'sig001a';


load([data_dir  task_file(1:8) '_3-spatial_analysis_results.mat']);
load([data_dir  task_file(1:8) '-Saccade_Direction_and_Amplitude_Analysis.mat']);
load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);


H = define_spatial_filter(filter_width);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

%---generic parameters---%
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
image_on_twin = 500;%don't use data within 500 ms of fixation onset due to strong visual response
min_window_width = 50; %time around peak to take for firing rate, based on FWHM of place cell fixation analysis
numshuffs = 10000; %number of times to shuffled/permutate saccade direction analysis
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
ITIstart_code = 15; %start of ITI/trial
img_on_code = 23; %image turned on
img_off_code = 24; %image turned off
smval = 30 ;%15 ms std, want to smooth at high frequency modulations
min_fix_dur = 100; %100 ms, don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva, don't want mini/micro saccades too small and hard to detect
min_saccades_with_spikes = 0.05;% 5% to remove neurons that have basically no activity since they pass too :(...
%essentially 0.4 Hz threshold if window is 50 ms in wide

%---saccade direction and amplitude parameters---%
bin_deg = 3;%1 %number of degrees per bin
bin_deg2 = 45; %initial degree bins to determine time of peak direction modualtion
smval_deg = 6;%18 %9 degrees std
bin_amplitude = 2; %dva for calculating saccade amplitude tuning
bin_amplitude2 = 4;%dva for estimating window of interest
max_amplitude = 16;%dva, approximately 95th percentile, don't have many large ones so remove

twinad1 = 200;%presaccade window
twinad2 = 400;%pos saccade window
median_sac_dur = 44;
tm = -twinad1+1:twinad2;
start_window = twinad1-min_fix_dur;%how early before saccade can you look for direction tuning
end_window = twinad1+min_fix_dur+median_sac_dur;%how late after saccade can you look for direction tuning,
%some neurons have large latencies


degrees = [0:bin_deg:360]-180;
degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Grab Saccade Aligned Activity--%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sac_aligned = saccade_aligned_firing{unit}; %fixation aligned firing rate
sac_dirs = saccade_directions{unit}; %saccade directions organized the same way as sac_algined
sac_amps = saccade_amplitudes{unit}/24; %saccade amplitudes organized the same way as sac_algined
fix_starts = fixation_starts{unit};
fix_ends = fixation_ends{unit};

%---Remove Counfounding Eye Movements for Spatially Modulated Neurons---%
%remove saccades in2in or out2out since in2out or out2in
%could be biased by field location creating artificial direction tuning
%but only do this if spatial (both criterion)
if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
    sac_aligned = sac_aligned(sac_in_out{unit} == 2 | sac_in_out{unit} == 4,:); %saccades in2in or out2out
    sac_dirs = sac_dirs(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for saccades in2in or out2out only
    sac_amps = sac_amps(sac_in_out{unit} == 2 | sac_in_out{unit} == 4); %directions for saccades in2in or out2out only
    fix_starts = fix_starts((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
    fix_ends = fix_ends((sac_in_out{unit} == 2 | sac_in_out{unit} == 4));
end

%also remove saccades that are too big since won't have many anyway
sac_aligned_amp = sac_aligned;
too_large = find(sac_amps > max_amplitude);
sac_aligned_amp(too_large,:) = [];
sac_amps(too_large)=[];
sac_fix_starts = fix_starts;
sac_fix_starts(too_large) = [];
sac_fix_ends = fix_ends;
sac_fix_ends(too_large) = [];

%---Modify Window if Necessary based on Fixation Durations and Remove Fixation Durations that are too Short---%
window_width = length(all_direction_windows{unit});
window_end = all_direction_windows{unit}(end);
window_start = all_direction_windows{unit}(1);
if window_end > end_window
    window_end = window_end-twinad1;
    median_fixation_dur = median(fix_ends);
    %if more than half the fixations are shorter than window length
    %cut window down to median fixation duration othwerise we will just
    %cut fixations shorter than window end
    if median_fixation_dur < window_end
        window = all_direction_windows{unit};
        window(window > median_fixation_dur+twinad1) = [];
        if length(window) < min_window_width;
            all_direction_windows{unit} = [];
        else
            all_direction_windows{unit} = window;
        end
        window_end = all_direction_windows{unit}(end)-twinad1;
    end
    %remove fixations shorter than window as these could be
    %contaminated by the next saccade
    fixations_too_short = find(fix_ends < window_end);
    sac_dirs(fixations_too_short) = [];
    sac_aligned(fixations_too_short,:) = [];
    fix_ends(fixations_too_short) = [];
    fix_starts(fixations_too_short) = [];
end
if window_start < twinad1-min_fix_dur
    fixations_too_short = find((twinad1+fix_starts) > window_start);
    sac_dirs(fixations_too_short) = [];
    sac_aligned(fixations_too_short,:) = [];
    fix_ends(fixations_too_short) = [];
    fix_starts(fixations_too_short) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Saccade Direction Plots---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
[~,si] = sort(sac_dirs);
sac_aligned_sorted = sac_aligned(si,:);
[trial,time] = find(sac_aligned_sorted == 1);
plot(time-twinad1,trial,'.k')
xlim([-twinad1 twinad2])
if ~isempty(trial)
    yl = [0 max(trial)+1];
    ylim([0 max(trial)+1]);
    hold on
    if ~isempty(all_direction_windows{unit})
        h = fill([all_direction_windows{unit}(1) all_direction_windows{unit}(end) ...
            all_direction_windows{unit}(end) all_direction_windows{unit}(1) all_direction_windows{unit}(1)]-twinad1,...
            [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
        %set(h,'facealpha',.25,'EdgeColor','None')
    end
    hold off
end
xlim([-twinad1 twinad2]);
ylabel('Ranked Sac. Dir.')
xlabel('Time from Saccade Start (ms)')
title('Fixation Aligned-Saccade Direction')
box off


%---Firing Rate Map---%
figure
firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
h = imagesc(firing_rate_map);
set(h,'alphadata',~isnan(firing_rate_map));
axis off
axis equal
colormap('jet')
colorbar
clim = caxis;
caxis([clim(1) maxfr])
title_str = ['\\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2)];
title(sprintf(title_str));

%---Firing Rate Curves---%
figure
[~,prefered_dirs,anti_prefered_dirs,smoothed_direction_curve] = ...
    select_prefred_indeces(direction_binned_firing_rate_curves{1,unit},degrees,sac_dirs,smval_deg,bin_deg);
hold on
dofill(tm,sac_aligned(prefered_dirs,:),'green',1,smval); %plot fixation aligned activity for prefered direction
dofill(tm,sac_aligned(anti_prefered_dirs,:),'black',1,smval); %plot fixation aligned activity for anti-prefered direction
yl = ylim;
if yl(1) < 0;
    yl(1) = 0;
    ylim([0 yl(2)]);
end
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-twinad1 twinad2]);
xlabel('Time from Saccade Start (ms)')
ylabel('Firing Rate (Hz)')
title('All Saccade Aligned Activity')


%---Plot Distribution of Saccade Directions and Preffered Direction---%
prefferedDirs = sort(sac_dirs(prefered_dirs));
antiprefferedDirs = sort(sac_dirs(anti_prefered_dirs));
if max(diff(prefferedDirs)) < 45
    prefLimits = [min(prefferedDirs) max(prefferedDirs)];
else
    disp('break in gaps...verify this works!')
    gaps = find(diff(prefferedDirs) > 45);
    if length(gaps) > 1
        error('wtf')
    end
    prefLimits = [prefferedDirs(1) prefferedDirs(gaps); prefferedDirs(gaps+1) prefferedDirs(end)];
end
if max(diff(antiprefferedDirs)) < 45
    anitLimits = [min(antiprefferedDirs) max(antiprefferedDirs)];
else
     disp('break in gaps...verify this works!')
    gaps = find(diff(antiprefferedDirs) > 45);
    if length(gaps) > 1
        error('wtf')
    end
    anitLimits = [antiprefferedDirs(1) antiprefferedDirs(gaps); antiprefferedDirs(gaps+1) antiprefferedDirs(end)];
end

[map, name, desc] = cmap('C6');
figure
subplot(1,2,1)
imagesc([0 1],sort(sac_dirs),sort(sac_dirs)');
hold on
for p = 1:size(prefLimits,1)
    plot([2 2],prefLimits(p,:),'g')
end
for a = 1:size(anitLimits,1)
    plot([2 2],anitLimits(a,:),'k')
end
hold off
xticks([]);
ylim([-180 180])
yticks([-180 -135 -90 -45 0 45 90 135 180])
box off
colormap(map)

subplot(1,2,2)
[counts,edges] = histcounts(sac_dirs,-180:8:180);
counts2 = [counts(end:-1:1) counts counts];
counts2 = filtfilt(1/smval_deg/2*ones(1,smval_deg/2),1,counts2);
counts2 = counts2(length(counts)+1:length(counts)*2);
plot(counts2',edges(1:end-1))
hold on
%plot(counts',edges(1:end-1)) %raw data
for p = 1:size(prefLimits,1)
    plot([0.1 0.1],prefLimits(p,:),'g')
end
for a = 1:size(anitLimits,1)
    plot([0.1 0.1],anitLimits(a,:),'k')
end
hold off
set(gca, 'YAxisLocation', 'Right','xdir','reverse')
ylim([-180 180])
yticks([-180 -135 -90 -45 0 45 90 135 180])
xticks([]);
box off

%for visualization inspection to make sure consistent with population
%level-data
figure
polar(edges/180*pi,[counts2 counts2(1)])


%---Polar Plots of Firing Rates---%
figure(101)
p = polarplot(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b');
title(sprintf(['Max FR. ' num2str(max(smoothed_direction_curve),2) ' Hz, MRL: ' num2str(mrls.all_saccades(unit),2) ' (' num2str(mrls.all_saccades_shuffled_prctile(unit),3) '%%)']))
rlim([0 max(smoothed_direction_curve)])


if (spatial_info.shuffled_rate_prctile(unit) > 95) && (spatial_info.spatialstability_halves_prctile(1,unit) > 95) && ~isnan(mrls.out2out(unit))
    %only run if spatially modulated or possibly modulated (95% skaggs OR 95% corr 1/2)
    %---Out2Out Fixations---%
    sac_out_out = saccade_aligned_firing{unit}(sac_in_out{unit} == 4,:); %fixation aligned firing rate
    dir_out_out = saccade_directions{unit}(sac_in_out{unit} == 4); %saccade directions
    fix_starts = fixation_starts{unit}(sac_in_out{unit} == 4);
    fix_ends = fixation_ends{unit}(sac_in_out{unit} == 4);
    
    %---Remove Fixations with Durations that are too short---%
    window_end = all_direction_windows{unit}(end);
    window_start = all_direction_windows{unit}(1);
    if window_end > end_window
        fixations_too_short = find(fix_ends < (window_end-twinad1));
        dir_out_out(fixations_too_short) = [];
        sac_out_out(fixations_too_short,:) = [];
        fix_ends(fixations_too_short) = [];
        fix_starts(fixations_too_short) = [];
    end
    if window_start < twinad1-min_fix_dur
        fixations_too_short = find((twinad1+fix_starts) > window_start);
        dir_out_out(fixations_too_short) = [];
        sac_out_out(fixations_too_short,:) = [];
        fix_ends(fixations_too_short) = [];
        fix_starts(fixations_too_short) = [];
    end
    
    %---Smoothed Plot of Firing Rate by Saccade Direction Out2Out ixations---%
    [~,~,~,smoothed_direction_curve] = ...
        select_prefred_indeces(direction_binned_firing_rate_curves{3,unit},degrees,dir_out_out,smval_deg,bin_deg);
    figure(101)
    hold on
    polarplot(degrees,[smoothed_direction_curve(end) smoothed_direction_curve],'b')
    hold off
    title(sprintf(['Max FR. ' num2str(max(smoothed_direction_curve),2) ' Hz, MRL: ' num2str(mrls.all_saccades(unit),2) ' (' num2str(mrls.all_saccades_shuffled_prctile(unit),3) '%%)' ...
        '\n Out2Out mrl: ' num2str(mrls.out2out(unit),2) ' (' num2str(mrls.out2out_shuffled_prctile(unit),3) '%%)']))
end


%% Sacade Amplitude Stuff
%
% %plots fixation aligned rasters so that can export to eps
%
% clar
% data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
% %data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Resorted\';
%
% %non-place cell amplitude tuning 1
% task_file = 'TO151208_3';
% unit_name = 'sig003b';
%
% bin_amplitude = 2; %dva for calculating saccade amplitude tuning
% bin_amplitude2 = 2;%dva for estimating window of interest
% max_amplitude = 16;%dva, approximately 95th percentile, don't have many large ones so remove
%
%
% load([data_dir task_file(1:8) '-Eyemovement_Locked_List_results.mat']);
% load([data_dir task_file(1:end) '-spatial_analysis_results.mat'],...
%     'spatial_info','eyepos','spike_times','binsize','Fs','filter_width')
% load([data_dir task_file(1:8) '-Saccade_amplitude_Analysis.mat']);
%
% H = define_spatial_filter(filter_width);
%
%
% unit_stats = unit_names;
% this_unit = [];
% for unit = 1:size(unit_stats,2)
%     if strcmpi(unit_stats{1,unit},unit_name)
%         this_unit = unit;
%         break
%     end
% end
%
% imageX = 800;
% imageY = 600;
%
%
% sac_amps = fixation_information{unit}(:,6); %fixation aligned firing rate
% sac_amps = sac_amps/24; %convert to dva
% fix_aligned = fixation_locked_firing{unit}; %saccade amplitudes
%
% %remove fixation within the first 500 ms of images on
% %amplitude is affected by time within image period
% time_from_image_on = fixation_information{unit}(:,4); %fixation aligned firing rate
% too_early = find(time_from_image_on < image_on_twin);
% fix_aligned(too_early,:) = [];
% sac_amps(too_early,:) = [];
%
% window_width = 100;
%
% %remove saccades that are too big since wont have many anyway
% too_large = find(sac_amps > max_amplitude);
% fix_aligned(too_large,:) = [];
% sac_amps(too_large)=[];
%
% fr = nandens(fix_aligned,smval,'gauss',1000,'nanflt'); %firing rate curve aligned to fixations
% window = all_windows{unit};
%
%
%  firing_rates = sum(fix_aligned(:,window),2)*1000/window_width;
%
%
% t = -twin1:twin2-1; %time for plotting
%
%
% %---Plot Firing Rate Curves for Large vs Small Amplitude Saccades---%
% small = prctile(sac_amps,25);
% large = prctile(sac_amps,75);
%
% figure
% hold on
% dofill(t,fix_aligned(sac_amps <= small,:),'black',1,smval); %smoothed fiirng rate curve
% dofill(t,fix_aligned(sac_amps >= large,:),'green',1,smval); %smoothed fiirng rate curve
% yl = ylim;
% if yl(1) < 0
%     ylim([0 yl(2)]);
%     yl(1) = 0;
% end
% plot([0 0],[yl(1) yl(2)],'k--')
% hold off
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% xlabel('Time From Fixation Start')
% ylabel('Firing Rate')
% legend('Small','Large')
% title(['Firing Rate Curves for Small and Large Saccades']);
% yl = ylim;
%
%
% figure
%
% %---Fixation Aligned Firing Rate Curve--%
% dofill(t,fix_aligned,'black',1,smval); %smoothed fiirng rate curve
% hold on
% h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],'r'); %window of interest
% % set(h,'facealpha',.25,'EdgeColor','None')
% plot([0 0],[yl(1) yl(2)],'k--')
% hold off
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% xlabel('Time from Fixation Start (ms)')
% ylabel('Firing Rate (Hz)')
% title('All Fixation Aligned Activity')
% ylim(yl)
%
%
% %---Fixation Aligned Raster Sorted by Saccade amplitude for out2out & in2in fixations---%
% figure
% [~,si] = sort(sac_amps);
% fix_aligned_sorted = fix_aligned(si,:);
% [trial,time] = find(fix_aligned_sorted == 1);
% plot(time-twin1,trial,'.k')
% xlim([-twin1 twin2])
% if ~isempty(trial)
%     ylim([0 max(trial)+1]);
%     hold on
%     h = fill([window(1) window(end) window(end) window(1) window(1)]-twin1,...
%         [0 0 max(trial)+1 max(trial)+1 0],'r');
%     set(h,'facealpha',.25,'EdgeColor','None')
%     hold off
% end
% xlim([-twin1 twin2]);
% set(gca,'Xtick',[-twin1 0 twin1 twin2])
% ylabel('Ranked Sac. Amp.')
% xlabel('Time from Fixation Start (ms)')
% title('Fixation Aligned-Saccade amplitude')
% box off
%
% %---Firing Rate Curve for Saccade Amplitude---%
% amps = [2:bin_amplitude:max_amplitude];
% figure
% plot(amps,binned_firing_rate_curves{unit})
% xlabel('Saccade Ampltiude')
% ylabel('Firing Rate')
% title(sprintf(['\\rho_{amp} = '  num2str(amplitude_correlations(unit),3) ' (' num2str(amplitude_correlations_percentile(unit),3) '%%)']))
% box off
%
% %---Firing Rate Map---%
% figure
% firing_rate_map = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all'); %get firing rate map
% maxfr = prctile(firing_rate_map(:),97.5); %set 97.5 percentile of firing rate map to help visualize
% h = imagesc(firing_rate_map);
% set(h,'alphadata',~isnan(firing_rate_map));
% axis off
% axis equal
% colormap('jet')
% colorbar
% clim = caxis;
% caxis([clim(1) maxfr])
% title_str = ['All images, peak rate = ' num2str(max(firing_rate_map(:)),3) ' Hz\n'];
% if spatial_info.shuffled_rate_prctile(unit) > 95;
%     title_str = [title_str 'Bits = ' num2str(spatial_info.rate(unit),3) ...
%         '(' num2str(spatial_info.shuffled_rate_prctile(unit),3) '%%) '];
% end
% if (spatial_info.spatialstability_halves_prctile(1,unit) > 95)
%     title_str = [title_str '\n \\rho_{1/2} = ' num2str(spatial_info.spatialstability_halves(1,unit),2) ...
%         '(' num2str(spatial_info.spatialstability_halves_prctile(1,unit),3) '%%)'];
% end
% if (spatial_info.spatialstability_even_odd_prctile(1,unit) > 95)
%     title_str = [title_str ' \\rho_{e/o} = ' num2str(spatial_info.spatialstability_even_odd(1,unit),2) ...
%         '(' num2str(spatial_info.spatialstability_even_odd_prctile(1,unit),3) '%%)'];
% end
% title(sprintf(title_str));

%% Visual & Memory Responses
% %plots fixation aligned rasters so that can export to eps

set(gcf,'renderer','Painters')

data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';

task_file = 'TO160114_3';
unit_name = 'sig004b';

load([data_dir  task_file(1:8) '-ListSQ-Visual_Response_Memory_results.mat']);
load([data_dir  task_file(1:8) '-ListSQ-Visual_Response_results.mat']);

this_unit = find(contains(unit_names,unit_name));
if isempty(this_unit)
    error('no unit found')
end


t12 = -twin1:twin2-1;
t13 = -twin1:twin3-1;
t133 = -twin1:twin3-1+250;
t3 = -twin3:twin3-1;
t14 = -twin1:twin4-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Visual Response---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Plot Firing rate locked to Image On 1 second window---%
figure
subplot(2,2,1)
hold on
[~,~,~,y3] =dofill(t12,time_lock_firing{this_unit,3},'black',1,smval); %smoothed image onset curve
plot([-twin1 twin2],[baseline_firing_rate(1,this_unit) baseline_firing_rate(1,this_unit)],'k--')%pre-image "baseline" line
hold off
ylabel('Firing Rate (Hz)')
xlabel('Time from Image On (ms)')
xlim([-twin1 twin2])
if epoch_data.rate_prctile(this_unit,3) > 90 || epoch_data.temporalstability_prctile(this_unit,3) > 90
    title(['bit ' num2str(epoch_data.rate_prctile(this_unit,3),3) '% ' ...
        '\rho = ' num2str(epoch_data.temporalstability(this_unit,3),2) ' ' ...
        num2str(epoch_data.temporalstability_prctile(this_unit,3),3)  '%' ])
end
box off

%---Plot raster for image on 1 second window---%
%in case we want to look at latency or something
subplot(2,2,3)
[trial,time] = find(time_lock_firing{this_unit,3} == 1); %novel
if ~isempty(trial)
    plot(time-twin1,trial,'.k')
    xlim([-twin1 twin2])
    ylim([0 max(trial)+1]);
end
ylabel('Trial #')
xlabel('Time from Image On (ms)')
box off

%---Firing rate locked to Image On 5 second window---%
subplot(2,2,2)
[~,~,~,y4] =dofill(t14,time_lock_firing{this_unit,5},'black',1,smval2); %smoothed firing rate curve
hold on
plot([-twin1 twin4],[baseline_firing_rate(1,this_unit) baseline_firing_rate(1,this_unit)],'k--')%pre-image "baseline" line
hold off
ylabel('Firing Rate (Hz)')
xlabel('Time from Image On (ms)')
xlim([-twin1 twin4])
if epoch_data.rate_prctile(this_unit,5) > 90 || epoch_data.temporalstability_prctile(this_unit,5) > 90
    title(['bit ' num2str(epoch_data.rate_prctile(this_unit,5),3) '% ' ...
        '\rho = ' num2str(epoch_data.temporalstability(this_unit,5),2) ' ' ...
        num2str(epoch_data.temporalstability_prctile(this_unit,5),3)  '%' ])
end
box off

%---Raster for Image On 5 second window---%
subplot(2,2,4)
[trial,time] = find(time_lock_firing{this_unit,5} == 1);
if ~isempty(trial)
    plot(time-twin1,trial,'.k')
    xlim([-twin1 twin4])
    ylim([0 max(trial)+1]);
end
ylabel('Trial #')
xlabel('Time from Image On (ms)')
box off

%---Scale plots to be the same and plot sig. difference from pre-image "baseline"---%
ymin = 0.85*min([y3 y4]);
if ymin < 0.1
    ymin = 0;
end
ymax = 1.2*max([y3 y4]);
if ymin ~= ymax
    
    %short 1 second window
    subplot(2,2,1)
    hold on
    ylim([ymin ymax])
    gaps = findgaps(find(( sig_visual_response(this_unit,1:twin1+twin2))));
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
    
    %long 5 second window
    subplot(2,2,2)
    hold on
    ylim([ymin ymax])
    gaps = findgaps(find((sig_visual_response(this_unit,:))));
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
end

n_str = [' Image Onset, n_ =' num2str(size(time_lock_firing{this_unit,1},1))];
if multiunit(this_unit)
    subtitle(['Multiunit ' unit_names{this_unit} n_str]);
else
    subtitle(['' unit_names{this_unit} n_str]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Memory Response---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Plot raster for image on-short window---%
subplot(2,3,1)
[trial,time] = find(time_lock_firing{this_unit,3}(nvr{this_unit} == 1 | nvr{this_unit}== 2,:) == 1); %all pairs
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
[trial,time] = find(time_lock_firing{this_unit,3}(nvr{this_unit} == 1,:) == 1); %novel
if ~isempty(trial)
    plot(time-twin1,trial,'.b')
    xlim([-twin1 twin2])
    ylim([0 max(trial)+1]);
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(time_lock_firing{this_unit,3}(nvr{this_unit} == 2,:) == 1); %repeat
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
dofill(t12,time_lock_firing{this_unit,3},'black',1,smval); %all trials
dofill(t12,time_lock_firing{this_unit,3}(nvr{this_unit} == 1,:),'blue',1,smval); %novel trials
dofill(t12,time_lock_firing{this_unit,3}(nvr{this_unit} == 2,:),'red',1,smval); %repeat trials
hold off
ylabel('Firing Rate (Hz)')
xlabel('Time from Image On (ms)')
xlim([-twin1 twin2])
if epoch_data.rate_prctile(this_unit,3) > 95 && epoch_data.temporalstability_prctile(this_unit,3) > 95
    title(['bit ' num2str(epoch_data.rate_prctile(this_unit,3),3) '% ' ...
        '\rho = ' num2str(epoch_data.temporalstability(this_unit,3),2) ' ' ...
        num2str(epoch_data.temporalstability_prctile(this_unit,3),3)  '%' ])
end
ylims(1,:) = ylim;
box off

%---Long 5 Second Window---%
subplot(2,3,4)
[trial,time] = find(time_lock_firing{this_unit,5}(nvr{this_unit} == 1 | nvr{this_unit}== 2,:) == 1);%all pairs
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
[trial,time] = find(time_lock_firing{this_unit,5}(nvr{this_unit} == 1,:) == 1); %novel
if ~isempty(trial)
    plot(time-twin1,trial,'.b')
    xlim([-twin1 twin4])
    ylim([0 max(trial)+1]);
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(time_lock_firing{this_unit,5}(nvr{this_unit} == 2,:) == 1); %repeat
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
dofill(t14,time_lock_firing{this_unit,5},'black',1,smval2); %all trials
dofill(t14,time_lock_firing{this_unit,5}(nvr{this_unit} == 1,:),'blue',1,smval2); %novel trials
dofill(t14,time_lock_firing{this_unit,5}(nvr{this_unit} == 2,:),'red',1,smval2); %repeat trials
hold off
ylabel('Firing Rate (Hz)')
xlabel('Time from Image On (ms)')
xlim([-twin1 twin4])
if epoch_data.rate_prctile(this_unit,5) > 95 && epoch_data.temporalstability_prctile(this_unit,5) > 95
    title(['bit ' num2str(epoch_data.rate_prctile(this_unit,5),3) '% ' ...
        '\rho = ' num2str(epoch_data.temporalstability(this_unit,5),2) ' ' ...
        num2str(epoch_data.temporalstability_prctile(this_unit,5),3)  '%' ])
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
gaps = findgaps(find(sig_short{this_unit}));
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
gaps = findgaps(find(sig_long{this_unit}));
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

subtitle(['Visual Response/Memory ' task_file(1:8) '_' unit_names{this_unit} ', n_ = ' num2str(sum(nvr{this_unit} == 1))]);