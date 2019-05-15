function [estimated_prefered_direction,prefered_dirs,anti_prefered_dirs,smoothed_firing] = select_prefred_indeces(binned_firing,binned_directions,dirs,smval,bin_size)
orginal_dirs = dirs; %store for later since may be modifying dirs variable
bf = [binned_firing(end-(3*smval):end) binned_firing binned_firing(1:3*smval)];%buffer so don't get edege artifacts from filtering
binned_firing = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
binned_firing = binned_firing(3*smval+2:end-3*smval);%remove buffers
smoothed_firing = binned_firing;%for output

%---Calculate Prefered Direction---%
%over smoothed so small peaks don't interfer detecting prefered direction
binned_firing2 = nandens(bf,3*smval,'gauss',1);
binned_firing2 = binned_firing2(3*smval+2:end-3*smval);%remove buffers
prefered_index = find((binned_firing2(1:end-1) == max(binned_firing2(1:end-1))));
prefered_index = prefered_index(1);
prefered_direction = 180/pi*binned_directions(prefered_index);%highest firing rate
estimated_prefered_direction = prefered_direction;

%---Calculate Anti-prefered Direction---%
%over smoothed so small peaks don't interfer detecting prefered direction
binned_firing2 = nandens(bf,3*smval,'gauss',1);
binned_firing2 = binned_firing2(3*smval+2:end-3*smval);%remove buffers
anti_index = find(binned_firing2(1:end-1) == min(binned_firing2(1:end-1)));
anti_index = anti_index(1);
anti_prefered_direction = 180/pi*binned_directions(anti_index); %lowest firing rate

%---Calculate FWHM for Prefered Direction---%
binned_firing2 = binned_firing-mean(binned_firing);
if abs(prefered_direction) > 90
    %rotate axis if prefered direction is near window edge
    zeroind = find(binned_directions == 0);
    binned_firing2 = [binned_firing2(zeroind:end) binned_firing2(1:zeroind-1)];
    
    binned_firing3 = nandens(bf,3*smval,'gauss',1);
    binned_firing3 = binned_firing3(3*smval+2:end-3*smval);%remove buffers
    binned_firing3 = [binned_firing3(zeroind:end) binned_firing3(1:zeroind-1)];
    prefered_index = find((binned_firing3(1:end-1) == max(binned_firing3(1:end-1))));
    prefered_index = prefered_index(1);
end

binned_firing2 = binned_firing2-mean(binned_firing2);
large_ind = find(binned_firing2 > binned_firing2(prefered_index)/2);%1/2 the max
gaps = findgaps(large_ind);
window = [];
for g =1:size(gaps,1)
    gp = gaps(g,:);
    gp(gp == 0) = [];
    if any(gp == prefered_index)
        window = gp;
        break
    end
end
deg_win = length(window)/2*bin_size;
if deg_win < 15
    deg_win = 15;
elseif deg_win > 45
    deg_win = 45;
end

%---Calculate FWHM for Anti-Prefered Direction---%
if abs(anti_prefered_direction) > 90
    %rotate axis if anti-prefered direction is near window edge
    zeroind = find(binned_directions == 0);
    binned_firing2 = [binned_firing2(zeroind:end) binned_firing2(1:zeroind-1)];
    
    binned_firing3 = nandens(bf,3*smval,'gauss',1);
    binned_firing3 = binned_firing3(3*smval+2:end-3*smval);%remove buffers
    binned_firing3 = [binned_firing3(zeroind:end) binned_firing3(1:zeroind-1)];
    anti_index = find(binned_firing3(1:end-1) == min(binned_firing3(1:end-1)));
    anti_index  = anti_index(1);
end

binned_firing2 = binned_firing2-mean(binned_firing2);
small_ind = find(binned_firing2 < binned_firing2(anti_index)/2);%1/2 the min
gaps = findgaps(small_ind);
window = [];
for g =1:size(gaps,1)
    gp = gaps(g,:);
    gp(gp == 0) = [];
    if any(gp == anti_index)
        window = gp;
        break
    end
end
deg_win_anti = length(window)/2*bin_size;
if deg_win_anti < 15
    deg_win_anti = 15;
elseif deg_win_anti > 45
    deg_win_anti = 45;
end

%---Select Directions in Prefered Direction---%
if prefered_direction > 180-deg_win %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif prefered_direction < -180+deg_win
    prefered_direction = prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
prefered_dirs = find((dirs <= prefered_direction+deg_win)& (dirs >= prefered_direction-deg_win));

%---Select Directions in Anti-Prefered Direction---%
dirs = orginal_dirs;%restore dirs variable since may have just modifid above
if anti_prefered_direction > 180-deg_win_anti %will have to look at negative angles too
    dirs(dirs < 0) = dirs(dirs < 0)+360;
elseif anti_prefered_direction < -180+deg_win_anti
    anti_prefered_direction = anti_prefered_direction+360;
    dirs(dirs < 0) = dirs(dirs < 0)+360;
else
    dirs = dirs;
end
anti_prefered_dirs = find((dirs <= anti_prefered_direction+deg_win_anti)& (dirs >= anti_prefered_direction-deg_win_anti));
end