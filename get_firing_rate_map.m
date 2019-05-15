function [ratemaps,timemaps,split_eyepos,split_spike_times] = get_firing_rate_map(trial_data,imageX,imageY,binsize,H,Fs,type)
%Written by Seth Koenig 9/19/2016 taken from various other codes.
%Code calculates firing rate map from eye position (dot position in cvtnew)
%and spike locations.
%
% Inputs:
%   1) trial_data:
%       a) trial_data{1}: position with rows by trial and columns by time
%             -even and odd rows are x and y position, respectively
%       b) trial_data{2}: spike times with rows by trial and columns by
%       time paralleling trial_data{1} in structure
%   2) imageX: horizontal size of display e.g. 800
%   3) imageY: vertical size of display e.g. 600
%   4) binsize: desired bin size of rate map e.g. 12 pixels
%   5) H: filter for smoothing rate map
%   6) Fs: sampling rate e.g. 1000 Hz
%   7) type: rate map type
%       a)'all': computes rate map for all trials
%       b)'half': computes rate maps for first and 2nd half of trials
%       c)'even_odd': computes rate map for even and odd trials
%
% Outputs:
%   1) ratemaps: firing rate maps
%   2) timemaps: coverage maps showing total time spent in a location
%   3) split_eyepos: eye position for halves or even/odds
%   4) split_spike_times: spike times for halves or even/odds
%   ***for half and even_odd outputs are cell arrays
%
% Code rechecked for bugs October 18, 2016 SDK

type = lower(type);
if strcmpi(type,'all')
    timemaps = filter_time(trial_data{1},imageX,imageY,Fs,binsize,H); %amount of time occupied over space
    filtered_space = filter_space(trial_data{1},trial_data{2},imageX,imageY,binsize,H); %spike count over space
    ratemaps = filtered_space./timemaps; %observed firing rate over space
    split_eyepos = [];
    split_spike_times = [];
    
elseif strcmpi(type,'half')
    ntrials = floor(size(trial_data{2},1)/2); %number of trials/2

    spike_times1 = trial_data{2}(1:ntrials,:); %spike times for first half
    spike_times2 = trial_data{2}(ntrials+1:ntrials*2,:);%spike times for second half
    
    eyepos1 = trial_data{1}(1:ntrials*2,:);%eye position for first half
    eyepos2 = trial_data{1}(ntrials*2+1:ntrials*2*2,:); %eye position for second half
    
    %rate map for first half
    timemaps1 = filter_time(eyepos1,imageX,imageY,Fs,binsize,H);
    filtered_space1 = filter_space(eyepos1,spike_times1,imageX,imageY,binsize,H);
    ratemap1 = filtered_space1./timemaps1;
    
    %rate map for second half
    timemaps2 = filter_time(eyepos2,imageX,imageY,Fs,binsize,H);
    filtered_space2 = filter_space(eyepos2,spike_times2,imageX,imageY,binsize,H);
    ratemap2 = filtered_space2./timemaps2;
    
    ratemaps{1} = ratemap1;
    ratemaps{2} = ratemap2;
    timemaps{1} = timemaps1;
    timemaps{2} = timemaps2;
    split_eyepos{1} = eyepos1;
    split_eyepos{2} = eyepos2;
    split_spike_times{1} = spike_times1;
    split_spike_times{2} = spike_times2;
    
elseif strcmpi(type,'even_odd')
    
    %eye position for even and odd trials
    x = trial_data{1}(1:2:end,:);
    y = trial_data{1}(2:2:end,:);
    
    x_odd = x(1:2:end,:);
    x_even = x(2:2:end,:);
    y_odd = y(1:2:end,:);
    y_even = y(2:2:end,:);
    
    %put back in to origianl format with x on odd rows and y on even rows
    eye_even = NaN(2*size(x_even,1),size(x_even,2));
    eye_even(1:2:end,:) = x_even;
    eye_even(2:2:end,:) = y_even;
    
    eye_odd = NaN(2*size(x_odd,1),size(x_odd,2));
    eye_odd(1:2:end,:) = x_odd;
    eye_odd(2:2:end,:) = y_odd;
    
    %spike times for even and odd trials
    spike_times_even = trial_data{2}(2:2:end,:);
    spike_times_odd = trial_data{2}(1:2:end,:);
    
    %rate map for even trials
    timemaps_even = filter_time(eye_even,imageX,imageY,Fs,binsize,H);
    filtered_space_even = filter_space(eye_even,...
        spike_times_even,imageX,imageY,binsize,H);
    ratemap_even = filtered_space_even./timemaps_even;
    
    %rate map for odd trials
    timemaps_odd = filter_time(eye_odd,imageX,imageY,Fs,binsize,H);
    filtered_space_odd = filter_space(eye_odd,...
        spike_times_odd,imageX,imageY,binsize,H);
    ratemap_odd = filtered_space_odd./timemaps_odd;
    
    ratemaps{1} = ratemap_even;
    ratemaps{2} = ratemap_odd;
    timemaps{1} = timemaps_even;
    timemaps{2} = timemaps_odd;
    split_eyepos{1} = eye_even;
    split_eyepos{2} = eye_odd;
    split_spike_times{1} = spike_times_even;
    split_spike_times{2} = spike_times_odd;
    
else
    error('Ratemap type not recognized')
end

    function [filtered_space] = filter_space(eyepos,spike_times,imageX,imageY,binsize,H)
        %caluclate total spikes over space
        [firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY);
        filtered_space = bin2(firing_location,binsize,binsize);
        filtered_space = imfilter(filtered_space,H,'replicate');
        filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
    end

    function [filtered_time] = filter_time(eyepos,imageX,imageY,Fs,binsize,H)
        %calculate the total time spent at any locations in binned pixels
        spatial_time = time_per_pixel(eyepos,imageX,imageY,Fs);
        filtered_time = bin2(spatial_time,binsize,binsize);
        zero_bins = (filtered_time == 0);
        filtered_time = imfilter(filtered_time,H,'replicate');
        filtered_time(zero_bins) = NaN;
        filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
    end

    function [firing_location] = pixel_spike_location(eyepos,spike_times,imageX,imageY)
        %written by Seth Konig August, 2014
        firing_location = zeros(imageY,imageX);
        
        [trial,time] = find(spike_times == 1);
        spikeind = sub2ind(size(spike_times),trial,time);
        
        x = eyepos(1:2:end,:);
        y = eyepos(2:2:end,:);
        x = stretch(x,imageX);
        y = stretch(y,imageY);
        
        xs = x(spikeind);
        ys = y(spikeind);
        [xs,ys] = remove_nans(xs,ys);
        
        %probably can't use logical indexing fast here since need to add spikes to pixels
        for i = 1:length(xs);
            firing_location(ys(i),xs(i)) =  firing_location(ys(i),xs(i)) +1;
        end
        
    end

    function [spatial_time] = time_per_pixel(eyepos,imageX,imageY,Fs)
        % written by Seth Konig August, 2014
        spatial_time = zeros(imageY,imageX);
        x = eyepos(1:2:end,:);
        y = eyepos(2:2:end,:);
        x = stretch2(x,imageX);
        y = stretch2(y,imageY);
        
        %probably can't use logical indexing fast here since need to add spikes to pixels
        xy_ind = sub2ind(size(spatial_time),y,x);
        for i = 1:length(xy_ind);
            spatial_time(xy_ind(i)) =  spatial_time(xy_ind(i)) +1;
        end
        spatial_time = spatial_time/Fs; %convert from ms to sec
    end


    function vector = stretch(vector,maximum) %turn into vector but keep nans
        vector = vector(1:end);
        vector(vector < 1) = 1;
        vector(vector > maximum) = maximum;
    end

    function vector = stretch2(vector,maximum) %turn into vector but remove nans
        vector = vector(1:end);
        vector(isnan(vector)) = [];
        vector(vector < 1) = 1;
        vector(vector > maximum) = maximum;
    end
end