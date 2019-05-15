function [ratemap,filtered_time] = outside_ratemap(eyepos,Fs,binsize,H,spike_times)
% modified 'get_firing_rate_map' Seth Konig 9.23.16
% code rechecked for bugs 1/22/2017 SDK


filtered_time = filter_time_outside_image(eyepos,Fs,binsize,H);
filtered_space = filter_space_outside(eyepos,spike_times,binsize,H);
ratemap = filtered_space./filtered_time;


    function [filtered_time] = filter_time_outside_image(eyepos,Fs,binsize,H)
        %calculate the total time spent at any locaitons outside the image
        spatial_time = time_per_pixel_outside(eyepos,Fs);
        filtered_time = bin2(spatial_time,binsize,binsize);
        zero_bins = (filtered_time == 0);
        filtered_time = imfilter(filtered_time,H,'replicate');
        filtered_time(zero_bins) = NaN;
        filtered_time = filtered_time(end:-1:1,:); %flip so rightside up
    end

    function [filtered_space] = filter_space_outside(eyepos,spike_times,binsize,H)
        %caluclate total spikes over space
        [firing_location] = pixel_spike_location_outside(eyepos,spike_times);
        filtered_space = bin2(firing_location,binsize,binsize);
        filtered_space = imfilter(filtered_space,H,'replicate');
        filtered_space = filtered_space(end:-1:1,:); %flip so rightside up
    end

    function [spatial_time] = time_per_pixel_outside(eyepos,Fs)
        % written by Seth Konig August, 2014
        
        spatial_time = zeros(900,1200);
        
        %eye position much more than these are saturated and or blinks
        min_y = -124;
        max_y = 724;
        min_x = -199;
        max_x = 1000;
        
        %get x and y
        x = eyepos(1:2:end,:);
        y = eyepos(2:2:end,:);
        
        %vectorize
        x = x(:);
        y = y(:);
        
        %remove nans
        x(isnan(x)) = [];
        y(isnan(y)) = [];
        
        %remove positions where x is outside reliable range
        y(x < min_x) = [];
        x(x < min_x) = [];
        y(x > max_x) = [];
        x(x > max_x) = [];
        
        %remove positions where y is outside reliable range
        x(y < min_y) = [];
        y(y < min_y) = [];
        x(y > max_y) = [];
        y(y > max_y) = [];
        
        %add so all values are positive
        x = round(x+200);
        y = round(y+150);
        
        %probably can't use logical indexing fast here since need to add spikes to pixels
        xy_ind = sub2ind(size(spatial_time),y,x);
        for i = 1:length(xy_ind);
            spatial_time(xy_ind(i)) =  spatial_time(xy_ind(i)) +1;
        end
        spatial_time = spatial_time/Fs; %convert from ms to sec
    end

    function [firing_location] = pixel_spike_location_outside(eyepos,spike_times)
        %written by Seth Konig August, 2014
        firing_location = zeros(900,1200);
        
        [trial,time] = find(spike_times == 1);
        spikeind = sub2ind(size(spike_times),trial,time);
        
        
        %eye position much more than these are saturated and or blinks
        min_y = -124;
        max_y = 724;
        min_x = -199;
        max_x = 1000;
        
        %get x and y
        x = eyepos(1:2:end,:);
        y = eyepos(2:2:end,:);
        
        %vectorize
        x = x(:);
        y = y(:);
        
        %remove positions where x is outside reliable range
        y(x < min_x) = NaN;
        x(x < min_x) = NaN;
        y(x > max_x) = NaN;
        x(x > max_x) = NaN;
        
        %remove positions where y is outside reliable range
        x(y < min_y) = NaN;
        y(y < min_y) = NaN;
        x(y > max_y) = NaN;
        y(y > max_y) = NaN;
        
        %add so all values are positive
        x = round(x+200);
        y = round(y+150);
        
        xs = x(spikeind);
        ys = y(spikeind);
        [xs,ys] = remove_nans(xs,ys);
        
        %probably can't use logical indexing fast here since need to add spikes to pixels
        for i = 1:length(xs);
            firing_location(ys(i),xs(i)) =  firing_location(ys(i),xs(i)) +1;
        end
        
    end
end