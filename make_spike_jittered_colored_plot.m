function make_spike_jittered_colored_plot(position,spike_times,src,subnum)
% make jittered spike plot but color codes by "time"
%src: subplot row and subplot col
%subnum: subplot number
%
% rechecked for bugs October 19, 2016 SDK

jitter = 12; %1/2 dva for spike locations

clrs = ['rb']; %gm want first vs second half

x = position(1:2:end,:);
y = position(2:2:end,:);

num_trials = size(spike_times,1);
step = floor(num_trials/length(clrs));

subplot(src(1),src(2),subnum)
plot(x',y','color',[0.8 0.8 0.8])
hold on
for group = 1:length(clrs)
    
    if group == 1
        trial_ind = 1:step;
    elseif group == length(clrs)
        trial_ind = (group-1)*step+1:num_trials;
    else
        trial_ind = (group-1)*step+1:group*step;
    end
    
    spk_times = spike_times(trial_ind,:);
    x_sub = x(trial_ind,:);
    y_sub = y(trial_ind,:);
    
    [trial,time] = find(spk_times == 1);
    spikeind = sub2ind(size(x_sub),trial,time);
    
    xs = x_sub(spikeind);
    ys = y_sub(spikeind);
    [xs,ys] = remove_nans(xs,ys);
    
    if ~isempty(xs)
        xs = xs+randi(jitter,length(xs),1);
        ys = ys+randi(jitter,length(ys),1);
        plot(xs,ys,['.' clrs(group)],'markersize',4)
    end
end
hold off
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Xcolor',[1 1 1]);
set(gca,'Ycolor',[1 1 1]);

xlim([-25 825])
ylim([-25 625])
axis equal
end
