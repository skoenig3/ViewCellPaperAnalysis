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

task = 'ListSQ';
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
Fs = 1000; %Hz sampling frequency
imageX = 800; %horizontal image size
imageY = 600; %vertical image size
image_on_twin = 500;

% time_locked_firing since ITI period is defined by 2 events (15 & 16)
image_on_code = 23;
image_off_code = 24;
reward_code = 3;
trial_start_code = 15; %ITI start

min_fix_dur = 100; %100 ms, don't want fixations that are too short since won't get a good idea of firing pattern
min_sac_amp = 48;%48 pixels = 2 dva, don't want mini/micro saccades too small and hard to detect

x_pos = [];
y_pos = [];
direction = [];
amplitude = [];

monkeys = {'Vivian','Tobii'};
for monk =2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir{1} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir{2} = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    for sess = 1:length(session_data)
        %read in task file data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{sess},task);
        if isempty(task_file)
            continue
        end
        
        disp(task_file(1:end-11))
        %load task file data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','hdr','fixationstats');
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %get important task specific information
        [itmlist,sequence_items,~] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [~,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,image_on_code);
        
        %preallocate space and parallel structure of cfg
        num_trials = length(cfg.trl); %total number of trials
        image_trials = zeros(1,num_trials); %image trials only
        for t = 1:num_trials %only take trials in which image was actually shown
            if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end);
                image_trials(t) = 1;
            end
        end
        
        
        %remove excess data associated with non image trials
        image_trials = find(image_trials);
        num_trials = length(image_trials); %total number of image trials
        fixationstats = fixationstats(image_trials); %only take the eye data from successful trials
        
        
        
        %1) trial #
        %2) %fixation start time from trial start
        %3) ordinal fixation #
        %4) fixation start relative to image onset
        %5) fixation duration
        %6) preceeding saccade amplitude
        %7) preceeding saccade duration
        %8) preceeding saccade direction
        %9) image novel/repeat
        
        xs = NaN(1,5000);
        ys = NaN(1,5000);
        dirs = NaN(1,5000);
        amps = NaN(1,5000);
        
        fix_ind = 1; %start fixation index @ 1
        for t = 1:num_trials %will use all trials since may be useful to look at eye movements later for all trials
            fixationtimes = fixationstats{t}.fixationtimes; %start and end times of fixations
            saccadetimes = fixationstats{t}.saccadetimes; % start and end times of saccades
            fixations = fixationstats{t}.fixations; %mean fixation posiiton
            xy = fixationstats{t}.XY;
            
            trial_start = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == trial_start_code); %time at start of trial
            imgon = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == image_on_code)-trial_start; %time of image on relative to start of trial
            imgoff = cfg.trl(image_trials(t)).alltim(cfg.trl(image_trials(t)).allval == image_off_code)-trial_start; %time of image off relative to start of trial
            
            % if monkey isn't paying attention data isn't probably
            % worth much plus have to cut off somewhere
            if imgoff-imgon > 1.5*imgdur-1
                imgoff = imgon+1.5*imgdur-1;
            end
            
            imgon = imgon+image_on_twin;%not analyzing data within first 500 ms of image presentation
            
            %find fixations and saccades that did not occur during the image period;
            %should also take care of the 1st fixation on the crosshair
            
            %fixation started before image turned on
            invalid= find(fixationtimes(1,:) < imgon);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            %fixation ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(fixationtimes(2,:) > imgoff);
            fixationtimes(:,invalid) = [];
            fixations(:,invalid) = [];
            
            
            %saccade started before image turned on
            invalid= find(saccadetimes(1,:) < imgon);
            saccadetimes(:,invalid) = [];
            
            %Saccaded ended after the image turned off so firing rate could corrupted by image turning off
            invalid= find(saccadetimes(2,:) > imgoff);
            saccadetimes(:,invalid) = [];
            
            for f = 1:size(fixationtimes,2);
                prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
                if isempty(prior_sac) %trial ended or eye outside of image
                    continue %try next one
                end
                sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
                if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough
                    
                    xs(fix_ind) = fixations(1,f);
                    ys(fix_ind) = fixations(2,f);
                    dirs(fix_ind) =  atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                        xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%saccade direction
                    amps(fix_ind) = sacamp;
                    
                    fix_ind = fix_ind+1;
                end
            end
        end
        
        xs(isnan(xs)) = [];
        ys(isnan(ys)) = [];
        dirs(isnan(dirs)) = [];
        amps(isnan(amps)) = [];
        
        x_pos = [x_pos xs];
        y_pos = [y_pos ys];
        direction = [direction dirs];
        amplitude = [amplitude amps];
        
    end
end

x_pos = round(x_pos);
x_pos(x_pos < 1) = 1;
x_pos(x_pos > imageX) = imageX;

y_pos = round(y_pos);
y_pos = imageY-y_pos; %flip
y_pos(y_pos < 1) = 1;
y_pos(y_pos > imageY) = imageY;

direction = direction*pi/180; %convert to radians
amplitude = amplitude/24;%convert to dva
%%


step = 50; %~2dva

x_bins = [0:step:imageX+step];
y_bins = [0:step:imageY+step];

observed_average_amplitude = NaN(length(y_bins)-2,length(x_bins)-2);
observed_std_amplitude = NaN(length(y_bins)-2,length(x_bins)-2);
observed_average_direction = NaN(length(y_bins)-2,length(x_bins)-2);
observed_mrl_direction = NaN(length(y_bins)-2,length(x_bins)-2);
for xi = 1:length(x_bins)-2
    for yi = 1:length(y_bins)-2
        ind = find((x_pos > x_bins(xi)) & (x_pos <= x_bins(xi+1)) & (y_pos > y_bins(yi)) & (y_pos <= y_bins(yi+1)));
        
        if ~isempty(ind)
            
            observed_average_amplitude(yi,xi) = mean(amplitude(ind));
            observed_std_amplitude(yi,xi) = std(amplitude(ind));
            
            observed_average_direction(yi,xi) = circ_mean(direction(ind));
            observed_mrl_direction(yi,xi) = circ_r(direction(ind));
            
        end
    end
end
%%

figure
subplot(2,2,1)
imagesc(observed_average_amplitude)
colormap('jet')
colorbar
title('Average Saccade Amplitude by Fixation Location')
axis off
axis equal

subplot(2,2,3)
imagesc(observed_std_amplitude)
colormap('jet')
colorbar
title('Standard Deviation of Saccade Amplitude by Fixation Location')
axis off
axis equal

subplot(2,2,2)
imagesc(observed_average_direction*180/pi)
colormap('jet')
colorbar
title('Average Saccade Direction by Fixation Location')
axis off
axis equal

subplot(2,2,4)
imagesc(1-observed_mrl_direction)
colormap('jet')
colorbar
title('Circular Variance (1-MRL) of Saccade Direction by Fixation Location')
axis off
axis equal

%%
[x,y] = meshgrid(1:size(observed_mrl_direction,2),1:size(observed_mrl_direction,1));

u = observed_mrl_direction.*cos(observed_average_direction);
v = -observed_mrl_direction.*sin(observed_average_direction);


figure
imagesc(observed_average_direction*180/pi)
hold on
quiver(x,y,u,v,'linewidth',1)
hold off
box off
axis off
%% Bootstrap to determine whether observed values are significant or not

numshuffs = 10000;

average_amplitude = NaN(length(y_bins)-2,length(x_bins)-2,numshuffs);
std_amplitude = NaN(length(y_bins)-2,length(x_bins)-2,numshuffs);
average_direction = NaN(length(y_bins)-2,length(x_bins)-2,numshuffs);
mrl_direction = NaN(length(y_bins)-2,length(x_bins)-2,numshuffs);

for shuff = 1:numshuffs
    
    shuffled_ind = randperm(length(x_pos));
    x_pos_shuff = x_pos(shuffled_ind);
    y_pos_shuff = y_pos(shuffled_ind);
    for xi = 1:length(x_bins)-2
        for yi = 1:length(y_bins)-2
            ind = find((x_pos_shuff > x_bins(xi)) & (x_pos_shuff <= x_bins(xi+1))...
                & (y_pos_shuff > y_bins(yi)) & (y_pos_shuff <= y_bins(yi+1)));
            
            if ~isempty(ind)
                
                average_amplitude(yi,xi,shuff) = mean(amplitude(ind));
                std_amplitude(yi,xi,shuff) = std(amplitude(ind));
                
                average_direction(yi,xi,shuff) = circ_mean(direction(ind));
                mrl_direction(yi,xi,shuff) = circ_r(direction(ind));
                
            end
        end
    end
end



sig_average_amplitude = NaN(length(y_bins)-2,length(x_bins)-2);
shuffled_average_amplitude = NaN(length(y_bins)-2,length(x_bins)-2);
sig_amplitude = NaN(length(y_bins)-2,length(x_bins)-2);
sig_average_direction = NaN(length(y_bins)-2,length(x_bins)-2);
sig_mrl_direction = NaN(length(y_bins)-2,length(x_bins)-2);
for xi = 1:length(x_bins)-2
    for yi = 1:length(y_bins)-2
        ind = find((x_pos > x_bins(xi)) & (x_pos <= x_bins(xi+1)) & (y_pos > y_bins(yi)) & (y_pos <= y_bins(yi+1)));
        
        if ~isempty(ind)
            
            %calculate observed percentile compared to shuffled
            shuffled_average_amplitude(yi,xi) = mean(squeeze(average_amplitude(yi,xi,:)));
            sig_average_amplitude(yi,xi) = 100*sum(observed_average_amplitude(yi,xi) > average_amplitude(yi,xi,:))/numshuffs;
            sig_std_amplitude(yi,xi) = 100*sum(observed_std_amplitude(yi,xi) > std_amplitude(yi,xi,:))/numshuffs;
            sig_mrl_direction(yi,xi) = 100*sum(observed_mrl_direction(yi,xi) > mrl_direction(yi,xi,:))/numshuffs;

            %percentiles don't work the same way for osberved direction
            sig_average_direction(yi,xi) = circ_mean(squeeze(average_direction(yi,xi,:)));
            
        end
    end
end

%%
figure
subplot(2,2,1)
imagesc(sig_average_amplitude)
colormap('jet')
colorbar
title('Bootstrapped Saccade Amplitude Percentile')
axis off
axis equal

subplot(2,2,3)
imagesc(sig_average_amplitude > 100-5/192)
colormap('jet')
colorbar
title('Bonferroni Corrected Sig/Not Sig')
axis off
axis equal

subplot(2,2,2)
imagesc(sig_mrl_direction)
colormap('jet')
colorbar
title('Bootstrapped Saccade Direction Percentile')
axis off
axis equal

subplot(2,2,4)
imagesc(sig_mrl_direction > 100-5/192)
colormap('jet')
colorbar
title('Bonferroni Corrected Sig/Not Sig')
axis off
axis equal



subitlte('Bootstrapped Data')

