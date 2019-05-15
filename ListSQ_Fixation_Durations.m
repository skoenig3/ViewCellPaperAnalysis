%written by Seth Konig September 14, 2016
%code extacts fixation durations, saccade rate, saccade_durations, and saccade amplitudes
%
% code rechecked 1/9/2017 SDK

clar %clear, clc
task = 'ListSQ';
ITIstart_code = 15; %start of ITI
img_on_code= 23; %start of image presentation
img_off_code = 24; %end of image presentation

fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
buffer = 100;

fltord = 60;
lowpasfrq = 100;
nyqfrq = 1000 ./ 2;
flt2 = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %100 Hz low pass filter

set = 0;%set index #
fixation_durations = cell(2,85);
saccade_amplitudes = cell(2,85);
saccade_durations = cell(2,85);
vel_profile = []; %algined to saccade start
smoothed_vel_profile = []; %smoothed algined to saccade start
smoothed_vel_profile2 = []; %smoothed more algined to saccade start
fixation_vel_profile = []; %algined to fixation start
fixation_smoothed_vel_profile = []; %smoothed algined to fixation start
fixation_smoothed_vel_profile2 = [];%smoothed more algined to fixation start
twin = 100;
which_monkey =[];
trial_count = [];


RMS_noise = cell(2,85); %variability in fixation position as a liberal estimate of eye tracking noise
%row 1 is x, row 2 is y
for monkey =1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 20]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for session = 1:length(session_data)
        task_file=get_task_data(session_data{session},task);
        if isempty(task_file)
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%
        %---Import Data---%
        %%%%%%%%%%%%%%%%%%%
        set = set +1; %set index+1
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'cfg','fixationstats','item_file','cnd_file');
        
        %get important task specific information
        [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);
        
        
        %set the image duration
        if str2double(task_file(3:8)) < 140805
            imgdur = 7000;
        else
            imgdur = 5000;
        end
        
        %---Prealoccate Memory---%
        fixation_durations{1,set} = NaN(96,40);
        fixation_durations{2,set} = NaN(96,40);
        saccade_amplitudes{1,set} = NaN(96,40);
        saccade_amplitudes{2,set} = NaN(96,40);
        saccade_durations{1,set} = NaN(96,40);
        saccade_durations{2,set} = NaN(96,40);
        RMS_noise{1,set} = NaN(192,40);%x
        RMS_noise{2,set} = NaN(192,40);%y
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---Get Behavioral Stats---%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_trials = length(cfg.trl); %number of trials
        
        %velocity aligned to saccade
        v_profile = NaN(4000,twin*2);
        sv_profile = NaN(4000,twin*2);
        sv_profile2 = NaN(4000,twin*2);
        
        %velocity aligned to fixation
        fv_profile = NaN(4000,twin*2);
        fsv_profile = NaN(4000,twin*2);
        fsv_profile2 = NaN(4000,twin*2);
        
        vind = 1;
        for t = 1:num_trials
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                
                %---get trial information---%
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code); %start time of trial
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turned on
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turned off
                
                % if monkey isn't paying attention data isn't probably
                % worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1
                    imgoff = imgon+1.5*imgdur-1;
                end
                
                %get novel/repeat information
                img_index = find(img_cnd == cfg.trl(t).cnd); %get image index
                if any(isnan(which_img(img_index)))%presentation error so skip trial, see get_image_numbers.m
                    continue
                end
                which_image = which_img(img_index); %image #
                nvr = novel_vs_repeat(img_index); %whether image was novel or repeat
                
                
                %---get fixation/saccaqde information---%
                fixationtimes = fixationstats{t}.fixationtimes; %fixtaion start and end times
                fixations = fixationstats{t}.fixations; %fixation locations
                saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                xy = fixationstats{t}.XY;
                
                % parse the input to remove the breaks in the eye data caused by looking
                % outside as indicated by the presence of NaNs
                vel = [];
                svel = [];
                svel2 = [];
                
                parsed_eyedat = preparse(xy);
                
                for p = 1:length(parsed_eyedat)
                    if any(~isnan(parsed_eyedat{p}(1,:)))
                        
                        %raw velocity
                        x = parsed_eyedat{p}(1,:);
                        y = parsed_eyedat{p}(2,:);
                        velx = diff(x);
                        vely = diff(x);
                        v = sqrt(velx.^2+vely.^2);
                        v(end+1) = v(end);
                        vel = [vel v];
                        
                        %low pass filtered velocity 30 Hz cutoff
                        x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                        y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                        xss = filtfilt(flt,1,x);
                        yss = filtfilt(flt,1,y);
                        xss = xss(101:end-101); %remove buffer after filtering
                        yss = yss(101:end-101); %remove buffer after filtering
                        x = x(101:end-101); %remove buffer after filtering
                        y = y(101:end-101); %remove buffer after filtering
                        
                        svelx = diff(xss);
                        svely = diff(yss);
                        sv = sqrt(svelx.^2+svely.^2);
                        sv(end+1) = sv(end);
                        svel = [svel sv];
                        
                        %low pass filtered velocity 100 Hz cutoff
                        x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                        y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                        xss = filtfilt(flt2,1,x);
                        yss = filtfilt(flt2,1,y);
                        xss = xss(101:end-101); %remove buffer after filtering
                        yss = yss(101:end-101); %remove buffer after filtering
                        x = x(101:end-101); %remove buffer after filtering
                        y = y(101:end-101); %remove buffer after filtering
                        
                        svelx = diff(xss);
                        svely = diff(yss);
                        sv2 = sqrt(svelx.^2+svely.^2);
                        sv2(end+1) = sv2(end);
                        
                        svel2 = [svel2 sv2];
                        
                    else
                        vel = [vel NaN(1,size(parsed_eyedat{p},2))];
                        svel = [svel NaN(1,size(parsed_eyedat{p},2))];
                        svel2 = [svel2 NaN(1,size(parsed_eyedat{p},2))];
                    end
                end
                
                %remove fixations that started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %remove fixations started ending after image turned off
                invalid= find(fixationtimes(2,:) > imgoff);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                invalid= find(saccadetimes(2,:) > imgoff);
                saccadetimes(:,invalid) = [];
                
                fixdurs = diff(fixationtimes)+1; %calculate fixation duration
                fixation_durations{nvr,set}(which_image,1:size(fixationtimes,2)) = fixdurs;
                
                sacdurs = diff(saccadetimes)+1; %calculate sacccade duration
                saccade_durations{nvr,set}(which_image,1:size(saccadetimes,2)) = sacdurs;
                
                for f = 1:size(fixationtimes,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                    
                    RMS_noise{1,set}(img_index,f) = std(xy(1,fixationtimes(1,f):fixationtimes(2,f)))/24; %std/rms in dva
                    RMS_noise{2,set}(img_index,f) = std(xy(2,fixationtimes(1,f):fixationtimes(2,f)))/24; %std/rms in dva
                    
                    if f == 1
                        sacamp = sqrt(sum((fixations(:,1)-[400;300]).^2));
                        saccade_amplitudes{nvr,set}(which_image,f) = sacamp;
                    else
                        prior_sac = find(saccadetimes(2,:) == fixationtimes(1,f)-1);%next fixation should start immediately after
                        if isempty(prior_sac) %no prior saccade so was proabbly looking off screen
                            continue;
                        end
                        sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                        
                        saccade_amplitudes{nvr,set}(which_image,f) = sacamp;
                        
                        
                        v_profile(vind,:) = vel(saccadetimes(1,prior_sac)-twin:saccadetimes(1,prior_sac)+twin-1);
                        sv_profile(vind,:) =svel(saccadetimes(1,prior_sac)-twin:saccadetimes(1,prior_sac)+twin-1);
                        sv_profile2(vind,:) =svel2(saccadetimes(1,prior_sac)-twin:saccadetimes(1,prior_sac)+twin-1);
                        
                        fv_profile(vind,:) = vel(fixationtimes(1,f)-twin:fixationtimes(1,f)+twin-1);
                        fsv_profile(vind,:) =svel(fixationtimes(1,f)-twin:fixationtimes(1,f)+twin-1);
                        fsv_profile2(vind,:) =svel2(fixationtimes(1,f)-twin:fixationtimes(1,f)+twin-1);
                        vind = vind+1;
                    end
                end
            end
            
        end
        which_monkey(set) = monkey;
        trial_count(set) = num_trials;
        
        vel_profile = [vel_profile; nanmedian(v_profile) ];
        smoothed_vel_profile = [smoothed_vel_profile; nanmedian(sv_profile)];
        smoothed_vel_profile2 = [smoothed_vel_profile2; nanmedian(sv_profile2)];
        
        
        fixation_vel_profile = [fixation_vel_profile; nanmedian(fv_profile) ];
        fixation_smoothed_vel_profile = [fixation_smoothed_vel_profile; nanmedian(fsv_profile)];
        fixation_smoothed_vel_profile2 = [fixation_smoothed_vel_profile2; nanmedian(fsv_profile2)];
    end
end

%remove excess empty matrices
saccade_durations (:,set+1:end) = [];
fixation_durations(:,set+1:end) = [];
saccade_amplitudes(:,set+1:end) = [];
%% Fixuation Durations by Ordinal Fixation #
all_fix_durs = []; %all fixaiton durations
nov_fix_durs = NaN(size(fixation_durations,2),20); %median by set novel fixation durations
rep_fix_durs = NaN(size(fixation_durations,2),20); %median by set repeat fixation durations
for set = 1:size(fixation_durations,2);
    nov_durs = fixation_durations{1,set}; %novel images
    rep_durs = fixation_durations{2,set}; %repeat images
    
    %remove fixations shorter than 100 ms in duration since removed from analysis
    nov_durs(nov_durs < 100) = NaN;
    rep_durs(rep_durs < 100) = NaN;
    
    nov_fix_durs(set,:) = nanmedian(nov_durs(:,1:20)); %novel
    rep_fix_durs(set,:) = nanmedian(rep_durs(:,1:20)); %repeat
    
    all_fix_durs = [all_fix_durs; nov_durs; rep_durs]; %all fixation durations
end

%---Stats Test---%
p_wilx = signrank(median(nov_fix_durs),median(rep_fix_durs)); %Wilcoxon signed rank test for zero median between novel and repeated images
%nonparametric, repeated measures (across multiple fixations), analysis
p_vals = [];
for f = 1:size(nov_fix_durs,2)
    [~,p_vals(f)] = ttest(nov_fix_durs(:,f),rep_fix_durs(:,f)); %paired ttest, not sure if really valid but easy to interpret
end

%---Plot Results---%
figure
hold on
plot(nanmedian(nov_fix_durs))
errorb(1:20,nanmedian(nov_fix_durs),nanstd(nov_fix_durs)./sqrt(sum(~isnan(nov_fix_durs))),'color','b')
plot(nanmedian(rep_fix_durs),'r')
errorb(1:20,nanmedian(rep_fix_durs),nanstd(rep_fix_durs)./sqrt(sum(~isnan(rep_fix_durs))),'color','r')
for f = 1:size(nov_fix_durs,2)
    if p_vals(f) < 0.05/size(nov_fix_durs,2) %Bonferroni correction
        plot(f,215,'k*')
    end
end
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
xlim([0 21])
ylim([160 220])
axis square
legend('Novel','Repeat')
title(['PW and TO n_{sessions} = ' num2str(size(fixation_durations,2)) ', p_{Wilcoxon} = ' num2str(p_wilx,3)])

%% Distribution of All Fixation Durations
durs = all_fix_durs;
durs(durs > 400) = [];
figure
hist(durs,301)
hold on
plot([round(nanmedian(durs)) round(nanmedian(durs))],[0 2000],'k--')
hold off
xlabel('Fixation Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Fixation Durations \n'...
    'mean = ' num2str(round(nanmedian(durs))) ' ms median = ' num2str(round(nanmedian(durs))) ' ms']))
axis square

%% Distribution of All Saccade Durations
%novel/differences really don't exist (effect ~1-2 ms)
all_sac_durs = [];
for set = 1:size(saccade_durations,2);
    all_sac_durs = [all_sac_durs; saccade_durations{1,set}; saccade_durations{2,set}];
end

durs =  all_sac_durs;
durs(durs > 80) = [];
figure
hist(durs,70)
hold on
plot([round(nanmedian(durs)) round(nanmedian(durs))],[0 14000],'k--')
hold off
xlim([10 80])
xlabel('saccade Duration (ms)')
ylabel('Count')
title(sprintf(['Distribution of Sacccade Durations \n'...
    'mean = ' num2str(round(nanmedian(durs))) ' ms median = ' num2str(round(nanmedian(durs))) ' ms']))
axis square
box off
%% Distribution of Saccade Ampltiudes

nov_sac_amps = NaN(size(saccade_amplitudes,2),20); %median by set novel saccade_amplitude
rep_sac_amps = NaN(size(saccade_amplitudes,2),20); %median by set repeat saccade_amplitude
all_amplitudes = [];
for set = 1:size(saccade_amplitudes,2);
    nov_sac = saccade_amplitudes{1,set}/24; %novel  images and convert from pixel to dva
    rep_sac = saccade_amplitudes{2,set}/24; %repeat images and convert from pixel to dva
    
    all_amplitudes = [all_amplitudes nov_sac rep_sac];
    
    
    %remove saccades smaller than 2 dva
    nov_sac(nov_sac < 2) = NaN;
    rep_sac(nov_sac < 2) = NaN;
    
    
    nov_sac_amps(set,:) = nanmedian(nov_sac(:,1:20)); %novel
    rep_sac_amps(set,:) = nanmedian(rep_sac(:,1:20)); %repeat
end


%---Stats Test---%
p_wilx = signrank(median(nov_sac_amps),median(rep_sac_amps)); %Wilcoxon signed rank test for zero median between novel and repeated images
%nonparametric, repeated measures (across multiple fixations), analysiss
p_vals = [];
for f = 1:size(nov_sac_amps,2)
    [~,p_vals(f)] = ttest(nov_sac_amps(:,f),rep_sac_amps(:,f)); %paired ttest, not sure if really valid but easy to interpret
end

%---Plot Results---%
figure
hold on
plot(nanmedian(nov_sac_amps))
errorb(1:20,nanmedian(nov_sac_amps),nanstd(nov_sac_amps)./sqrt(sum(~isnan(nov_sac_amps))),'color','b')
plot(nanmedian(rep_sac_amps),'r')
errorb(1:20,nanmedian(rep_sac_amps),nanstd(rep_sac_amps)./sqrt(sum(~isnan(rep_sac_amps))),'color','r')
for f = 1:size(nov_sac_amps,2)
    if p_vals(f) < 0.05/size(nov_sac_amps,2) %Bonferroni correction
        plot(f,8.5,'k*')
    end
end
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
xlim([0 21])
ylim([5 9])
axis square
legend('Novel','Repeat')
title(['PW and TO n_{sessions} = ' num2str(size(saccade_amplitudes,2)) ', p_{Wilcoxon} = ' num2str(p_wilx,3)])
%% Post-hoc Power Analysis How many pairs of images you would need to show a difference

n = sampsizepwr('t2',[mean(nanmedian(nov_durs(:,3:12)')) std(nanmedian(nov_durs(:,3:12)'))],...
    [mean(nanmedian(rep_durs(:,3:12)')) std(nanmedian(rep_durs(:,3:12)'))],0.9)
%%

figure


%---Vivian---%
mns = mean(nov_fix_durs(which_monkey == 1,:),2);
[r,p] = corrcoef(1:length(mns),mns);
P = polyfit(1:length(mns),mns',1);
yfit = polyval(P,1:length(mns));

subplot(2,3,1)
hold on
plot(mean(nov_fix_durs(which_monkey == 1,:),1))
plot(mean(rep_fix_durs(which_monkey == 1,:),1))
hold off
legend('Novel','Repeat')
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('Vivian')
yl = ylim;

subplot(2,3,2)
plot(mns,'.k')
hold on
plot(1:length(mns),yfit,'k')
hold off
xlabel('Session #')
ylabel('Avg. Novel Fix. Dur. (ms)')
ylim(yl)
xlim([0 length(mns)+1])
box off
title(['Vivian: r = ' num2str(r(2),2) ', p = ' num2str(p(2),2) ', m = ' num2str(P(1),3)])

[r,p] = corrcoef(trial_count(which_monkey == 1),mns);
subplot(2,3,3)
plot(mns,trial_count(which_monkey == 1),'k.')
ylabel('# of Trials Completed')
xlabel('Avg. Novel Fix. Dur. (ms)')
title(['Vivian: r = ' num2str(r(2),2) ', p = ' num2str(p(2),2)])
box off

[r,p] = corrcoef(1:length(mns),trial_count(which_monkey == 1));
subplot(2,3,6)
plot(1:length(trial_count(which_monkey == 1)),trial_count(which_monkey == 1),'k.')
xlabel('Session #')
ylabel('Trial Count')
xlim([0 length(mns)+1])
box off
title(['Vivian: r = ' num2str(r(2),2) ', p = ' num2str(p(2),2)])


%---Tobii---%
mns = mean(nov_fix_durs(which_monkey == 2,:),2);
[r,p] = corrcoef(1:length(mns),mns);
P = polyfit(1:length(mns),mns',1);
yfit = polyval(P,1:length(mns));

subplot(2,3,4)
hold on
plot(mean(nov_fix_durs(which_monkey == 2,:),1))
plot(mean(rep_fix_durs(which_monkey == 2,:),1))
hold off
legend('Novel','Repeat')
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('Tobii')
yl = ylim;

subplot(2,3,5)
plot(mns,'.k')
hold on
plot(1:length(mns),yfit,'k')
hold off
xlabel('Session #')
ylabel('Avg. Novel Fix. Dur. (ms)')
ylim(yl)
xlim([0 length(mns)+1])
box off
title(['Tobii: r = ' num2str(r(2),2) ', p = ' num2str(p(2),2) ', m = ' num2str(P(1),3)])

subtitle('Fixation Durations by Session # for ListSQ Recording Sessions')
%%
figure
hold on
plot(-twin:twin-1,1000/24*(mean(vel_profile)))
plot(-twin:twin-1,1000/24*(mean(smoothed_vel_profile)))
plot(-twin:twin-1,1000/24*(mean(smoothed_vel_profile2)))

yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Velocity (dva/s)')
legend('Raw','LPF 30 Hz','LPF 100 Hz')

%%
raw = mean(vel_profile);
raw = raw-mean(raw(1:twin));
raw = raw/max(raw);

filtered = mean(smoothed_vel_profile);
filtered = filtered-mean(filtered(1:twin));
filtered = filtered/max(filtered);

filtered2 = mean(smoothed_vel_profile2);
filtered2 = filtered2-mean(filtered2(1:twin));
filtered2 = filtered2/max(filtered2);

figure
hold on
plot(-twin:twin-1,raw)
plot(-twin:twin-1,filtered);
plot(-twin:twin-1,filtered2);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlabel('Time from Saccade Start (ms)')
ylabel('Normalized Velocity')
legend('Raw','LPF 30 Hz','LPF 100 Hz')

%%
all_x_noise = [];
all_y_noise = [];
set_mean_noise = [];

for set = 1:length(RMS_noise)
    if ~isempty(RMS_noise{1,set})
        all_x_noise = [all_x_noise; RMS_noise{1,set}(:)];
        all_y_noise = [all_y_noise; RMS_noise{2,set}(:)];
        
        set_mean_noise = [set_mean_noise [nanmedian(RMS_noise{1,set}(:)); nanmedian(RMS_noise{2,set}(:))]];
    end
end
%%
x_99 = prctile(all_x_noise,99);
y_99 = prctile(all_y_noise,99);
%%
%%
figure
subplot(1,2,1)
hold on
histogram(set_mean_noise(1,:),25,'EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor','b')
histogram(set_mean_noise(2,:),25,'EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor','r')
hold off
legend('X','Y')
xlabel('Average Noise (\sigma)')
ylabel('Session Count')
title('Session Average')

subplot(1,2,2)
hold on
histogram(all_x_noise,100,'EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor','b')
histogram(all_y_noise,100,'EdgeAlpha',0.5,'FaceAlpha',0.5,'FaceColor','r')
hold off
xlabel('Individual Fixation Noise (\sigma)')
ylabel('Fixation Count')
xlim([0 1])
title(['99 pertcentile X = ' num2str(x_99,3) ', 99 percentile Y = ' num2str(y_99,3)])
%%
tm = -twin:twin-1;
figure
hold on
plot(tm,24*nanmedian(fixation_vel_profile))
plot(tm,24*nanmedian(vel_profile))
hold off
xlabel('Time from Eye Movement Start (ms)')
ylabel('Eye Velocity (dva/sec)')
xlim([-75 75])
title('Raw Eye velocity')

figure
hold on
plot(tm,24*nanmedian(fixation_smoothed_vel_profile))
plot(tm,24*nanmedian(smoothed_vel_profile))
hold off
xlabel('Time from Eye Movement Start (ms)')
ylabel('Eye Velocity (dva/sec)')
xlim([-75 75])
title('Over Smoothed')


figure
hold on
plot(tm,24*nanmedian(fixation_smoothed_vel_profile2))
plot(tm,24*nanmedian(smoothed_vel_profile2))
hold off
xlabel('Time from Eye Movement Start (ms)')
ylabel('Eye Velocity (dva/sec)')
xlim([-75 75])
title('Less Smoothed')


%% Seperate Fixation Durations by "high" and "low" recognition

%---Based on Change in Fixation Duration---%
%under idea that large changes are correlated with high recogntion trials
nov_small_diff = [];%low
rep_small_diff = [];%low
nov_large_diff = [];%high
rep_large_diff = [];%high

%---Based soley on Repeat Fixation Duration---%
%under idea that long repeat fixation durations are correlated with high recogntion trials
nov_short_rep = [];%low
rep_short_rep = [];%low
nov_long_rep = [];%high
rep_long_rep = [];%high

%---Based soley on Novel Fixation Duration---%
%under idea that long repeat fixation durations are correlated with high recogntion trials
nov_short_nov = [];%low
rep_short_nov = [];%low
nov_long_nov = [];%high
rep_long_nov = [];%high

for set = 1:size(fixation_durations,2);
    nov_durs = nanmedian(fixation_durations{1,set}(:,1:18)'); %novel images (:,3:12)
    rep_durs = nanmedian(fixation_durations{2,set}(:,1:18)'); %repeat images (:,3:12)

    %---Based on Change in Fixation Duration---%
    diff = rep_durs-nov_durs;
    
    upper_33 = find(diff >= prctile(diff,67));%high
    lower_33 = find(diff <= prctile(diff,33));%low
    
    nov_small_diff = [nov_small_diff; nanmedian(fixation_durations{1,set}(lower_33,1:18))]; %low
    rep_small_diff = [rep_small_diff; nanmedian(fixation_durations{2,set}(lower_33,1:18))];%low
    nov_large_diff = [nov_large_diff; nanmedian(fixation_durations{1,set}(upper_33,1:18))];%high
    rep_large_diff = [rep_large_diff; nanmedian(fixation_durations{2,set}(upper_33,1:18))];%high

    
    %---Based soley on Repeat Fixation Duration---%
    upper_33 = find(rep_durs >= prctile(rep_durs,67));%high
    lower_33 = find(rep_durs <= prctile(rep_durs,33));%low
    
    nov_short_rep = [nov_short_rep; nanmedian(fixation_durations{1,set}(lower_33,1:18))]; %low
    rep_short_rep = [rep_short_rep; nanmedian(fixation_durations{2,set}(lower_33,1:18))];%low
    nov_long_rep = [nov_long_rep; nanmedian(fixation_durations{1,set}(upper_33,1:18))];%high
    rep_long_rep = [rep_long_rep; nanmedian(fixation_durations{2,set}(upper_33,1:18))];%high
    
    %---Based soley on Novel Fixation Duration---%
    upper_33 = find(nov_durs >= prctile(nov_durs,67));%high
    lower_33 = find(nov_durs <= prctile(nov_durs,33));%low
    
    nov_short_nov = [nov_short_nov; nanmedian(fixation_durations{1,set}(lower_33,1:18))]; %low
    rep_short_nov = [rep_short_nov; nanmedian(fixation_durations{2,set}(lower_33,1:18))];%low
    nov_long_nov = [nov_long_nov; nanmedian(fixation_durations{1,set}(upper_33,1:18))];%high
    rep_long_nov = [rep_long_nov; nanmedian(fixation_durations{2,set}(upper_33,1:18))];%high
end

figure
subplot(2,2,1)
hold on
plot(nanmedian(nov_small_diff),'b')
errorb(1:18,nanmedian(nov_small_diff),nanstd(nov_small_diff)./sqrt(sum(~isnan(nov_small_diff))),'color','b')
plot(nanmedian(rep_small_diff),'r')
errorb(1:18,nanmedian(rep_small_diff),nanstd(rep_small_diff)./sqrt(sum(~isnan(rep_small_diff))),'color','r')
plot(nanmedian(nov_large_diff),'g')
errorb(1:18,nanmedian(nov_large_diff),nanstd(nov_large_diff)./sqrt(sum(~isnan(nov_large_diff))),'color','g')
plot(nanmedian(rep_large_diff),'k')
errorb(1:18,nanmedian(rep_large_diff),nanstd(rep_large_diff)./sqrt(sum(~isnan(rep_large_diff))),'color','k')
hold off
xlim([0 20])
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Low Nov','Low Rep','High Nov','High Rep')
title('Seperated by Change in Fixation Duration')
ylim([150 250])

subplot(2,2,2)
hold on
plot(nanmedian(nov_short_rep),'b')
errorb(1:18,nanmedian(nov_short_rep),nanstd(nov_short_rep)./sqrt(sum(~isnan(nov_short_rep))),'color','b')
plot(nanmedian(rep_short_rep),'r')
errorb(1:18,nanmedian(rep_short_rep),nanstd(rep_short_rep)./sqrt(sum(~isnan(rep_short_rep))),'color','r')
plot(nanmedian(nov_long_rep),'g')
errorb(1:18,nanmedian(nov_long_rep),nanstd(nov_long_rep)./sqrt(sum(~isnan(nov_long_rep))),'color','g')
plot(nanmedian(rep_long_rep),'k')
errorb(1:18,nanmedian(rep_long_rep),nanstd(rep_long_rep)./sqrt(sum(~isnan(rep_long_rep))),'color','k')
hold off
ylim([150 250])
xlim([0 20])
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('Seperated by Repeat Fixation Duration')

subplot(2,2,3)
hold on
plot(nanmedian(nov_short_nov),'b')
errorb(1:18,nanmedian(nov_short_nov),nanstd(nov_short_nov)./sqrt(sum(~isnan(nov_short_nov))),'color','b')
plot(nanmedian(rep_short_nov),'r')
errorb(1:18,nanmedian(rep_short_nov),nanstd(rep_short_nov)./sqrt(sum(~isnan(rep_short_nov))),'color','r')
plot(nanmedian(nov_long_nov),'g')
errorb(1:18,nanmedian(nov_long_nov),nanstd(nov_long_nov)./sqrt(sum(~isnan(nov_long_nov))),'color','g')
plot(nanmedian(rep_long_nov),'k')
errorb(1:18,nanmedian(rep_long_nov),nanstd(rep_long_nov)./sqrt(sum(~isnan(rep_long_nov))),'color','k')
hold off
ylim([150 250])
xlim([0 20])
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('Seperated by Novel Fixation Duration')



 