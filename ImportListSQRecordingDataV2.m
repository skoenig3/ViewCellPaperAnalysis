function ImportListSQRecordingDataV2(data_dir,figure_dir,session_data)
% Version 2.0 of ImportListSQRecordingData. It was re-written slighlty to
% handle importing sesion data directly from an excel sheet. Seth Konig 1/7/2016
%
%  Code imports cch25f.sav data to calibrate eye data. Eye data and
%  LFP/single unit data is also imported from the .nex file. Finally code
%  implemetns Cluster Fix to detect fixation and saccades.
%
% Inputs:
%   1) data_dir: directory containing .nex files for session
%   2) session_data: data containing info about the session data
%   3) figure_dir: folder for saving figures
% Outputs:
%   1) saves preprocessed data to data_dir
%
% Code rechecked for bugs August, 2016

figure_dir = [figure_dir '\Calibration\'];

%grab important file information from session_data
task = 'ListSQ';
cch25_file = session_data.calibration_file;%cch25f.sav
[task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,waveform_count]...
    = get_task_data(session_data,task);
if isempty(task_file)
    warning('No ListSQ file could be found. Exiting function...')
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Import Plexon_nex File Data---%%%
% define the trials solely based on the marker events
cfg=[];
cfg.channel = {};
cfg.dataset       = [data_dir task_file];
cfg.trialfun      = 'trialfunListSQ';
cfg.datatype      = 'continuous';

warning('off','all') %supress isrow warning
cfg = ft_definetrial(cfg);
warning('on','all')

% read the data from file and preprocess them
% Want spike data (valid units only), LFPs, and X and Y eye data
fieldtripdefs
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
num_units = 0;
waveform_channels = [];
for l = 1:length(hdr.label);
    if ~isempty(strfind(hdr.label{l},'sig')); %all recorded "units"
        if isempty(strfind(hdr.label{l},'wf')); %all recorded "units" waveforms, dont want for now
            if isempty(strfind(hdr.label{l}(end),'i'));%i is for invalidated waveforms?, rest are valid waveforms
                cfg.channel = [cfg.channel hdr.label(l)];
                num_units = num_units+1;
            end
        else
            if isempty(strfind(hdr.label{l}(7),'i'));%i is for invalidated waveforms, rest are valid waveforms
                waveform_channels= [waveform_channels l];
            end
        end
    elseif length(hdr.label{l}) == 1 && ~isempty(strfind(hdr.label{l},'X')); %horiziontal eye data
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif length(hdr.label{l}) == 1 && ~isempty(strfind(hdr.label{l},'Y')); %vertical eye data
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(hdr.label{l},'pupil')) %pupil data does not exist for all recordings
        cfg.channel = [cfg.channel hdr.label(l)];
    elseif ~isempty(strfind(hdr.label{l},'AD'));  %LFP data by channel
        cfg.channel = [cfg.channel hdr.label(l)];
    end
end
data = getPlexonTrialData(cfg);

% import waveforms seperately for some reason couldn't get it to work above
waveforms = cell(2,length(waveform_channels));
dataWF=readNexFile(cfg.dataset);
valid_units = [];
for wv = 1:length(dataWF.neurons);
    if ~strcmpi(dataWF.neurons{wv}.name(end),'i') %i for invalid waveforms
        valid_units = [valid_units wv];
    end
end

for wv = 1:length(valid_units)
    waveforms{1,wv} = dataWF.waves{valid_units(wv)}.waveforms;%waveorm shapes
    waveforms{2,wv} = dataWF.waves{valid_units(wv)}.timestamps;%waveform timestamps
    
    unit = strfind(unit_names,dataWF.neurons{valid_units(wv)}.name);
    unit = find(~cellfun(@isempty,unit));
    if length(waveforms{2,wv}) ~= waveform_count(unit)
        emailme([task_file ' spike count does not match excel file'])

        error(['Warning number of waveforms imported differnet than expected: ' ...
            'Imported ' num2str(length(waveforms{2,wv})) ' but expected ' num2str(waveform_count(wv))]);
    end
end

%when files are cut then merged they can have gaps in the data so want to
%remove these gaps. They appear to be due to time stamp start/ends that
%didn't get fixed properly by plexUtil. SDK 1/10/16.
% Gaps can only be in continuous data since they're replaced with NaNs
% while spike channels are just buffed with 0's
lfp_channels = find_desired_channels(cfg,'LFP');
gaps = cell(1,13);
for channel = lfp_channels(1):length(data)
    for trl = 1:length(data(1).values);
        d = data(channel).values{trl};
        if any(isnan(d))
            gaps{channel} = [gaps{channel} [trl; sum(isnan(d))]];
        end
    end
end

gapsQ = cellfun(@isempty,gaps);
lfp_channels = find_desired_channels(cfg,'LFP');
if any(~gapsQ) %then there is a gap in the data somewhere
    warning('Gaps have been found in the data')
    reply = input('Do you wish to continue Y/N [Y]:','s');
    if strfind(lower(reply),'y')
        %first check that all the gaps are only on the continuous channels
        gap_channels = find(~gapsQ);
        if any(gap_channels < min(lfp_channels))
            error('Something is not right with the gaps in the data')
        end
        
        %second check that all the gaps are only on the continuous channels
        same_gaps = zeros(length(gap_channels));
        for  g = 1:length(gap_channels);
            for gg = 1:length(gap_channels);
                if gaps{gap_channels(g)} == gaps{gap_channels(g)}
                    same_gaps(g,gg) = 1;
                end
            end
        end
        
        if all(all(same_gaps)) %then all the gaps are the same which is good
            %then lets fix these trials
            gap_trials = gaps{gap_channels(1)}(1,:);
            for trial = 1:length(gap_trials)
                nan_ind = find(isnan(data(gap_channels(1)).values{gap_trials(trial)}));
                disp(['Fixing gaps in trial # ' num2str(gap_trials(trial)) '. Fixing ' num2str(length(nan_ind)) ' gaps!'])
                %first fix the data
                for chan = 1:length(data)
                    data(chan).values{gap_trials(trial)}(nan_ind) = []; %just remove them
                end
                %second fix the cfg just for those trials
                %adjust end time and timing of events. Could adjust other
                %trials but that would pointless as it would mean nothing
                %with this analysis
                cfg.trl(gap_trials(trial)).endsmpind = cfg.trl(trial).endsmpind-length(nan_ind);
                tim = cfg.trl(gap_trials(trial)).alltim;
                after = find(tim-cfg.trl(gap_trials(trial)).begsmpind > nan_ind(1));
                tim(after) = tim(after) - length(nan_ind);
                cfg.trl(gap_trials(trial)).alltim = tim;
            end
        else
            error('Something is not right with the gaps in the data')
        end
        
        disp('Gaps in data is now fixed')
    else
        error('Gaps not fixed session will now abort')
    end
end

%---Verify that trial data is the same length across events and spikes/LFPs
%may happen on the last trial(s) if ended recording in middle of trial should
%be only off by 1 trial but can occassionally happen on more for unknown reasons 
if length(cfg.trl) > length(data(1).values)
    disp('Warning: there is less data trials than trials defined by strobes')
    disp([num2str(length(cfg.trl)) ' trials in cfg while ' num2str(length(data(1).values)) ' in spike data'])
    reply = input('Do you wish to remove cfg.trl trials Y/N [Y]:','s');
    if strfind(lower(reply),'y')
        num_trials = length(data(1).values);
        cfg.trl = cfg.trl(1:num_trials);
    else
        error('Data structure and cfg structure do not match in length')
    end
    
elseif length(data(1).values) > length(cfg.trl)
    error('There are more data trials than trials defined by strobes')
end

disp('Data imported Successfully')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---calibrate eye data---%%%
imageX = 800;
imageY = 600;

[tform] = Calibrate_Plexon_EyeData(data_dir,cch25_file,figure_dir); %get the calibration function

eyechans = find_desired_channels(cfg,'eye');
outside_xy = cell(1,length(data(eyechans(1)).values)); %values outside image
% replaces outliers (more than 1 dva outside image with NaNs not sure
for t = 1:length(data(eyechans(1)).values);
    x = data(eyechans(1)).values{t}; %second to last index in structure array is horizontal eye data;
    y = data(eyechans(2)).values{t}; %last index in structure array is vertical eye data;
    outside_xy{t} = NaN(2,length(x));
    
    [x,y] = tformfwd(tform,x,y); %calibrate: transform votlages into "dva", 24 pixels/dva
    
    x = 24*x; %convert from cortex dva to pixels
    y = 24*y; %convert from cortex dva to pixels
    x = x+imageX/2;
    y = y+imageY/2;
    
    %get x-values outside image border
    outside_xy{t}(1,x < 1) = x(x < 1);
    outside_xy{t}(2,x < 1) = y(x < 1);
    outside_xy{t}(1,x > imageX) = x(x > imageX);
    outside_xy{t}(2,x > imageX) = y(x > imageX);
    
    y(x < 1) = NaN;
    x(x < 1) = NaN;
    y(x > imageX) = NaN;
    x(x > imageX)= NaN;
    
    %get y-values outsie image border
    outside_xy{t}(1,y < 1) = x(y < 1);
    outside_xy{t}(2,y < 1) = y(y < 1);
    outside_xy{t}(1,y > imageY) = x(y > imageY);
    outside_xy{t}(2,y > imageY) = y(y > imageY);
    
    x(y < 1) = NaN;
    y(y < 1) = NaN;
    x(y > imageY) = NaN;
    y(y > imageY) = NaN;
    
    %store calibrated data
    data(eyechans(1)).values{t} = x;
    data(eyechans(2)).values{t} = y;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Run Cluster Fix---%%%
% use Cluster Fix to detect fixations and saccades in XY eye data
% since timing is very important NaNs must remain. Trials with NaNs or with eye
% data less lasting less than 200 ms (less than 1 saccade and fixation) are ignored
% and trials are processed in chunks of valid eye data.
disp('Running Cluster Fix')
fixationstats = cell(1,length(data(end-1).values));
for t = 1:length(data(eyechans(1)).values);
    if ~isempty(data(end-1).values{t})
        fixationstats{t} = ClusterFix_Plexon([data(eyechans(1)).values{t};data(eyechans(2)).values{t}]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save the Data---%%%
save([data_dir task_file(1:end-11) '-preprocessed.mat'],'cch25_file','cfg','data',...
    'fixationstats','hdr','item_file','cnd_file','task_file','num_units','waveforms','tform',...
    'multiunit','unit_confidence','sorting_quality','outside_xy');
disp(['Data  from ' task_file ' Successfully Preprocessed and Saved'])
end