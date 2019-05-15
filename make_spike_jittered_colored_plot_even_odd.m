function make_spike_jittered_colored_plot_even_odd(position,spike_times,src,subnum)
% make jittered spike plot but color codes by "time"
%src: subplot row and subplot col
%subnum: subplot number
%
% rechecked for bugs October 19, 2016 SDK

jitter = 12; %1/2 dva for spike locations

clrs = ['rb']; %gm want first vs second half

x = position(1:2:end,:);
y = position(2:2:end,:);

x_odd = x(1:2:end,:);
y_odd = y(1:2:end,:);
spike_times_odd = spike_times(1:2:end,:);

y_even = y(2:2:end,:);
x_even = x(2:2:end,:);
spike_times_even = spike_times(2:2:end,:);

subplot(src(1),src(2),subnum)
plot(x',y','color',[0.8 0.8 0.8])
hold on

[trial,time] = find(spike_times_odd  == 1);
spikeind = sub2ind(size(x_odd),trial,time);

xs = x_odd(spikeind);
ys = y_odd(spikeind);
[xs,ys] = remove_nans(xs,ys);

if ~isempty(xs)
    xs = xs+randi(jitter,length(xs),1);
    ys = ys+randi(jitter,length(ys),1);
    plot(xs,ys,['.' clrs(1)],'markersize',4)
end

[trial,time] = find(spike_times_even  == 1);
spikeind = sub2ind(size(x_even),trial,time);

xs = x_even(spikeind);
ys = y_even(spikeind);
[xs,ys] = remove_nans(xs,ys);

if ~isempty(xs)
    xs = xs+randi(jitter,length(xs),1);
    ys = ys+randi(jitter,length(ys),1);
    plot(xs,ys,['.' clrs(2)],'markersize',4)
end

hold off
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Xcolor',[1 1 1]);
set(gca,'Ycolor',[1 1 1]);

xlim([-25 825])
ylim([-25 625])
axis equal
