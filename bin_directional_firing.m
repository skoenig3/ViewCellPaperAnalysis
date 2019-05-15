function [mean_binned_firing_rate,degrees,mrl] = bin_directional_firing(bin_deg,saccade_aligned_firing,time_window,saccade_directions)
%bin data into small degree bins
degrees = [0:bin_deg:360]-180;
binned_firing_rate = cell(1,length(degrees));
for bins = 2:length(degrees)
    these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
    binned_firing_rate{bins} = 1000*nansum(saccade_aligned_firing(these_dirs,time_window),2)/length(time_window);
end
degrees = degrees(2:end);
degrees = degrees*pi/180;
degrees = [degrees degrees(1)];
mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
mrl = circ_r(saccade_directions*pi/180,sum(saccade_aligned_firing(:,time_window),2)); %MRL for unbinned data

%NOT rechecked for bugs since currently not used
%binned data
% degrees = [0:3:360]-180; %binned degrees
% binned_firing_rate = cell(1,length(degrees));
% for bins = 2:length(degrees)
%     these_dirs = (saccade_directions < degrees(bins) & saccade_directions >= degrees(bins-1));
%     binned_firing_rate{bins} = 1000*nansum(saccade_aligned_firing(these_dirs,time_window),2)/length(time_window);
% end
% degrees = degrees(2:end);
% degrees = degrees*pi/180;
% degrees = [degrees degrees(1)];
% mean_binned_firing_rate = cellfun(@nanmean,binned_firing_rate(2:end));
% mbfr = mean_binned_firing_rate;
% mbfr(isnan(mbfr)) = 0;
%mrl = circ_r(degrees(1:end-1),mbfr,bin_deg/180*pi); %calculate for binned data
end
