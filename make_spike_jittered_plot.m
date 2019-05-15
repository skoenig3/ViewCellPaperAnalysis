function make_spike_jittered_plot(position,spike_times,src,subnum)
%src: subplot row and subplot col
%subnum: subplot number
%
% rechecked for bugs October 19, 2016 SDK

jitter = 12; %1/2 dva for spike locations

x = position(1:2:end,:);
y = position(2:2:end,:);

[trial,time] = find(spike_times == 1);
spikeind = sub2ind(size(x),trial,time);

xs = x(spikeind);
ys = y(spikeind);
[xs,ys] = remove_nans(xs,ys);

xs = xs+randi(jitter,length(xs),1);
ys = ys+randi(jitter,length(ys),1);

subplot(src(1),src(2),subnum)
hold on
plot(x',y','color',[0.8 0.8 0.8])
plot(xs,ys,'.r','markersize',4)
hold off
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Xcolor',[1 1 1]);
set(gca,'Ycolor',[1 1 1]);

xlim([-25 825])
ylim([-25 625])
axis equal
end