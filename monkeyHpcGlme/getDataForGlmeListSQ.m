function getDataForGlmeListSQ(monkey,session,data_dir,session_data)
%Function grabs data so that we can run GLME later.
%Function runs 2 types of anlaysis: 1) GLM on whole data, 2) GLME on
%fixation-aligned data.
%
%written by Seth Koeni 7/24/22

%---Processing and Task Parameters---%
task = 'ListSQ';
colors ='rgbm';
shapes = 'xo';

imageX = 800;
imageY = 600;
img_on_code = 23; %cortex code when image turns on
img_off_code = 24; %cortex code when image turns off
ITIstart_code = 15; %start of ITI/trial

min_fix_dur = 100; %100 ms %don't want fixations that are too short since won't get a good idea of firing pattern
max_fix_dur = 500; %500 ms, anything longer is probably more than 1 eye movement
min_sac_amp = 48;%48 pixels = 2 dva don't want mini/micro saccades too small and hard to detect

Fs = 1000; %Hz sampling frequency
fixwin = 5;%size of fixation window on each crosshair
smval = 30;%2*std of gaussian kernel so 15 ms standard deviation
min_blks = 2; %only analyzes units with at least 2 novel/repeat blocks (any block/parts of blocks)
twin1 = 200;%200 ms before fixation
twin2 = 400;%400 ms after start of fixation
image_on_twin = 500;%how much time to ignore eye movements for to remove strong visual response though some may last longer

%bin sizes for raw data
eyeTimeBinSize = 40;%in ms/samples since at 1000 Hz
numSacBins = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---import task and unit data---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task_file = get_task_data(session_data,task);
if isempty(task_file)
    warning('No file could be found for specificed task. Exiting function...')
    return
end
try
    load([data_dir task_file(1:end-11) '-spatial_analysis_results.mat']);
catch
    disp('Spatial Analysis File not found continueing')
    return
end
load([data_dir task_file(1:end-11) '-preprocessed.mat']);

%Save as new variableso can reload later...kind of dumb but that was how it was written
absolute_fixationstats = fixationstats;
absolute_cfg = cfg;

[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]...
    = get_task_data(session_data,task);
[multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
    multiunit,unit_confidence,sorting_quality);
clear unit_names

if num_units == 0 %if no units exit function
    disp([task_file(1:8) ': no units could be found. Exiting function...'])
    return;
end

valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
if all(isnan(valid_trials(1,:)))
    disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
    return %since no valid units skip analysis
end

%remove units with too few trials
%these are the absolute minimum data required to do data analysis
%set the image duration
if str2double(task_file(3:8)) < 140805
    imgdur = 7000;
else
    imgdur = 5000;
end

%get important task specific information
[itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
[which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code);

%filter parameters
H = define_spatial_filter(filter_width);


%---Pre-allocate space for Fixations during Sequence Trials Analysis---%
allFiringRateMaps = cell(1,num_units);
sigSpatialScores = NaN(1,num_units);
spatialCorr = NaN(1,num_units);
allFixationData = cell(1,num_units); %smoothed
allFixationLockedFiring = cell(1,num_units); %fixation-average firing rate, parallels allFixationData
allFixationLockedFiringRaw = cell(1,num_units); %unsmoothed, binary, does not match the size of other variables
allFixationLockedFiringBinned = cell(1,num_units); %raw 25 ms bins
for unit = 1:num_units
    %check if we have data
    if isempty(eyepos{unit})
        continue %no data for this neuron
    end
    
    %determine if significant spatial unit+
    spatialCorr(unit) = spatial_info.spatialstability_halves(unit);
    if (spatial_info.shuffled_rate_prctile(unit) > 95) ... %skagg 95%+
            && (spatial_info.spatialstability_halves_prctile(unit) > 95) %spatial stability
        sigSpatialScores(unit) = 1;
    else
        sigSpatialScores(unit) = 0;
    end
    
    %preallocate space to matrix and table
    fix_ind = 1; %fixation # so can track in variables above
    allFixationLockedFiring{unit} = NaN(3000,twin1+twin2 + 1); %smoothed spike trains locked to fixations
    allFixationLockedFiringRaw{unit} = NaN(3000,twin1+twin2 + 1); %raw spike trains locked to fixations
    allFixationData{unit} = table('Size',[3000,15],...
        'VariableTypes',{'double','double','cell',...
        'double','double','double',...
        'double','double','double','double','double',...
        'double','double',...
        'double','double'},...
        'VariableNames',{'monkey','session','sessionName',...
        'imageNumber','imageID','fixationNumber',...
        'meanResponse','posX','posY','saccadeDir','saccadeAmp',...
        'fixationTime','logFixationTimage',...
        'fixationDuration','novelRepeatImage'});
    
    %25 ms binneed data so differnet size than the rest
    fixIndBinned = 1;
    allFixationLockedFiringBinned{unit} = table('Size',[20e3,16],...
        'VariableTypes',{'double','double','cell',...
        'double','double','double','double',...
        'double','double','double','double','double',...
        'double','double',...
        'double','double'},...
        'VariableNames',{'monkey','session','sessionName',...
        'imageNumber','imageID','fixationNumber','eyePhase',...
        'meanResponse','posX','posY','saccadeDir','saccadeAmp',...
        'fixationTime','logFixationTimage',...
        'fixationDuration','novelRepeatImage'});
    
    
    
    %%%---Calculate Firing Rate Map---%%%
    allFiringRateMaps{unit} = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize*2,H,Fs,'all'); %get firing rate map
    
    
    %%%---Calculate Firing Rate Locked to Fixations---%%%
    fixationstats = absolute_fixationstats; %reload because written over below
    cfg = absolute_cfg; %reload because written over below
    num_trials = length(cfg.trl);%number of trials
    for t = 1:num_trials
        if t >= valid_trials(1,unit) && t <= valid_trials(2,unit) %only valid trials for this unit
            if any(cfg.trl(t).allval == img_on_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == ITIstart_code);
                imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-trial_start; %when image turns on
                imgoff = cfg.trl(t).alltim(cfg.trl(t).allval == img_off_code)-trial_start; %when image turns off
                
                % if monkey isn't paying attention and looked away image presentation
                % is now longer than imgdur (because of cumulative looking time)
                % so data isn't probably worth much plus have to cut off somewhere
                if imgoff-imgon > 1.5*imgdur-1 %cut off trial at 1.5x length of desired looking time
                    imgoff = imgon+1.5*imgdur-1;
                end
                imgon = imgon+image_on_twin; %to avoid visual response and strong central bias
                
                fixationtimes = fixationstats{t}.fixationtimes; %fixation start and end times
                saccadetimes = fixationstats{t}.saccadetimes; %saccade start and end times
                fixations = round(fixationstats{t}.fixations); %mean fixation location
                xy = fixationstats{t}.XY; %xy eye trace
                
                %find fiations and saccades that did not occur during the image period;
                %should also take care of the 1st fixation on the crosshair
                
                %fixation started before image turned on
                invalid= find(fixationtimes(1,:) < imgon);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %fixation started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(fixationtimes(1,:) > imgoff-twin2);
                fixationtimes(:,invalid) = [];
                fixations(:,invalid) = [];
                
                %saccade started before image turned on
                invalid= find(saccadetimes(1,:) < imgon);
                saccadetimes(:,invalid) = [];
                
                %saccade started after the image turned off and/or firing rate could corrupted by image turning off
                invalid= find(saccadetimes(1,:) > imgoff-twin2);
                saccadetimes(:,invalid) = [];
                
                %which image monkey viewed if nan or empty then skip cuz bad trial
                img_index = find(img_cnd == cfg.trl(t).cnd);
                if isempty(img_index) || any(isnan(which_img(img_index)))
                    continue
                end
                thisNovelRepeat = novel_vs_repeat(img_index);
                thisImageNumber = which_img(img_index);
                
                %get fixation locked response
                rawSpikes = data(unit).values{t};
                [~,spikes] = nandens(data(unit).values{t},smval,'gauss',Fs); %was find(data(unit).values{t}); %spike trains for this trial
                for f = 2:size(fixations,2)%ignore first fixation not sure where it was/possibly contaminated anyway
                    prior_sac = find(fixationtimes(1,f) == saccadetimes(2,:)+1);%next fixation should start immediately after
                    if isempty(prior_sac) %trial ended or eye outside of image
                        continue %try next one
                    end
                    sacamp = sqrt(sum((xy(:,saccadetimes(2,prior_sac))-xy(:,saccadetimes(1,prior_sac))).^2)); %saccade amplitude
                    fix_dur = fixationtimes(2,f)-fixationtimes(1,f)+1;%this fixation duration
                    if sacamp >= min_sac_amp && fix_dur >= min_fix_dur && fix_dur <= max_fix_dur %next fixation has to be long enough & Fixation large enough
                        fixt = fixationtimes(1,f)-imgon;
                        sacdir = atan2d(xy(2,saccadetimes(2,prior_sac))-xy(2,saccadetimes(1,prior_sac)),...
                            xy(1,saccadetimes(2,prior_sac))-xy(1,saccadetimes(1,prior_sac)));%saccade direction
                        
                        %store everything in the table
                        allFixationData{unit}.monkey(fix_ind) = monkey;
                        allFixationData{unit}.session(fix_ind) = session;
                        allFixationData{unit}.sessionName{fix_ind} = task_file(1:8);
                        allFixationData{unit}.imageNumber(fix_ind) = img_index;
                        allFixationData{unit}.imageID(fix_ind) = thisImageNumber;
                        allFixationData{unit}.fixationNumber(fix_ind) = f;
                        allFixationData{unit}.meanResponse(fix_ind) = mean(spikes(saccadetimes(1,prior_sac):fixationtimes(2,f)));%mean response fixation + prior saccade
                        allFixationData{unit}.posX(fix_ind) = fixations(1,f);
                        allFixationData{unit}.posY(fix_ind) = fixations(2,f);
                        allFixationData{unit}.saccadeDir(fix_ind) = sacdir;
                        allFixationData{unit}.saccadeAmp(fix_ind) = sacamp;
                        allFixationData{unit}.fixationTime(fix_ind) = fixt;
                        allFixationData{unit}.logFixationTimage(fix_ind) = log(fixt);
                        allFixationData{unit}.fixationDuration(fix_ind) = fix_dur;
                        allFixationData{unit}.novelRepeatImage(fix_ind) = thisNovelRepeat;
                        
                        %get fixation locked data
                        allFixationLockedFiring{unit}(fix_ind,:) = spikes(fixationtimes(1,f)-twin1:fixationtimes(1,f)+twin2);
                        allFixationLockedFiringRaw{unit}(fix_ind,:) = rawSpikes(fixationtimes(1,f)-twin1:fixationtimes(1,f)+twin2);
                        fix_ind = fix_ind+1; %add index +1
                        
                        
                        %binned firing rate representation
                        %assume saccades ~50 ms in duration
                        numBins = numSacBins + floor(fix_dur/eyeTimeBinSize);%for the saccade may be slight overlap here and there
                        startTime = fixationtimes(1,f)- numSacBins*eyeTimeBinSize;
                        binTimes = startTime:eyeTimeBinSize:startTime+eyeTimeBinSize*numBins;
                        for nB = 1:numBins-1
                            binResponseCount = sum(rawSpikes(binTimes(nB):binTimes(nB+1)-1));
                            allFixationLockedFiringBinned{unit}.monkey(fixIndBinned) = monkey;
                            allFixationLockedFiringBinned{unit}.session(fixIndBinned) = session;
                            allFixationLockedFiringBinned{unit}.sessionName{fixIndBinned} = task_file(1:8);
                            allFixationLockedFiringBinned{unit}.imageNumber(fixIndBinned) = img_index;
                            allFixationLockedFiringBinned{unit}.imageID(fixIndBinned) = thisImageNumber;
                            allFixationLockedFiringBinned{unit}.fixationNumber(fixIndBinned) = f;
                            allFixationLockedFiringBinned{unit}.eyePhase(fixIndBinned) = nB;
                            allFixationLockedFiringBinned{unit}.meanResponse(fixIndBinned) = binResponseCount;
                            allFixationLockedFiringBinned{unit}.posX(fixIndBinned) = fixations(1,f);
                            allFixationLockedFiringBinned{unit}.posY(fixIndBinned) = fixations(2,f);
                            allFixationLockedFiringBinned{unit}.saccadeDir(fixIndBinned) = sacdir;
                            allFixationLockedFiringBinned{unit}.saccadeAmp(fixIndBinned) = sacamp;
                            allFixationLockedFiringBinned{unit}.fixationTime(fixIndBinned) = fixt;
                            allFixationLockedFiringBinned{unit}.logFixationTimage(fixIndBinned) = log(fixt);
                            allFixationLockedFiringBinned{unit}.fixationDuration(fixIndBinned) = fix_dur;
                            allFixationLockedFiringBinned{unit}.novelRepeatImage(fixIndBinned) = thisNovelRepeat;
                            fixIndBinned = fixIndBinned + 1;
                        end
                    end
                end
            end
        end
    end
    
    %clean up
    allFixationLockedFiring{unit}(fix_ind:end,:) = [];
    allFixationLockedFiringRaw{unit}(fix_ind:end,:) = [];
    allFixationData{unit}(fix_ind:end,:) = [];
    allFixationLockedFiringBinned{unit}(fixIndBinned:end,:) = [];
end


%---Save the Data---%
save([data_dir  task_file(1:8) '_glmeData.mat'],'unit_stats','eyeTimeBinSize',...
    'allFiringRateMaps','spatialCorr','sigSpatialScores',...
    'allFixationData','allFixationLockedFiring','allFixationLockedFiringRaw',...
    'allFixationLockedFiringBinned');

end