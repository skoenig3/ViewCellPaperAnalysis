function make_spike_jittered_plot_threshold(position,spike_times,threshold_matrix,src,subnum)
% make jittered spike plot threhsolded spikes only
%src: subplot row and subplot col
%subnum: subplot number

jitter = 12; %1/2 dva

x = position(1:2:end,:);
y = position(2:2:end,:);

[trial,time] = find(spike_times == 1);
spikeind = sub2ind(size(x),trial,time);

xs = x(spikeind);
ys = y(spikeind);
[xs,ys] = remove_nans(xs,ys);

%for speed up need to put these on the edge if they're their otherwise they
%wont' fit with the matrix
xs(xs < 1) = 1;
xs(xs > 800) = 800;
ys(ys < 1) = 1;
ys(ys > 600)= 600;

xys = sub2ind(size(threshold_matrix),ys,xs);

threshold_matrix = threshold_matrix(end:-1:1,:);
xs(threshold_matrix(xys) == 0) = [];
ys(threshold_matrix(xys) == 0) = [];

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