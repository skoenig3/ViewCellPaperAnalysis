% Plots figures for publication by plotting things individuall

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\TO Recording Files\';

task_file = 'TO151208';
unit_name = 'sig003b';

load([data_dir  task_file(1:8) '-Place_Cell_Analysis.mat']);
load([data_dir  task_file(1:end-11) '-spatial_analysis_results.mat']);

this_unit = [];
for unit = 1:size(unit_stats,2)
    if strcmpi(unit_stats{1,unit},unit_name)
        this_unit = unit;
        break
    end
end

fix_locked_firing = list_fixation_locked_firing{this_unit};
fix_in_out = in_out{this_unit};
t = -twin1:twin2-1;

num_in = sum(fix_in_out == 1);
num_out = sum(fix_in_out == 4);
downsample = floor(num_out/num_in)/2; %many more fixations outside so downsample to show ~equal number
figure

%% Plot Firing rates and rasters for in vs out
%---Fixations in->out vs out->out---%

out_matrix = fix_locked_firing(fix_in_out == 4,:);
out_matrix = out_matrix(1:downsample:end,:);
[trial,time] = find(out_matrix == 1);
plot(time-twin1,(trial),'.b')
hold on
if ~isempty(trial)
    b4 = max(trial);
else
    b4 = 0;
end
[trial,time] = find(fix_locked_firing(fix_in_out == 1,:) == 1);
trial = trial+b4;
plot(time-twin1,(trial),'.r')
if ~isempty(trial)
    ylim([0 max(trial)])
else
    ylim([0 b4])
end
box off
plot([0 0],[0 max(trial)+1],'k--')
ylabel('Occurence')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
xlabel('Time from Fixation Start (ms)')
title('Fixation Aligned Rasters')
axis square

figure
hold on
[~,~,~,y_list,~] = dofill(t,fix_locked_firing(fix_in_out == 1,:),'red',1,smval);%out-> in
dofill(t,fix_locked_firing(fix_in_out == 4,:),'blue',1,smval);%out->out
[pks,locs] = findpeaks(y_list,'MinPeakWidth',40);
locs(pks < 0.66*max(y_list)) = [];
pks(pks < 0.66*max(y_list)) = [];
plot(locs-twin1,pks,'*k')
yl = ylim;
if yl(1) < 0
    yl(1) = 0;
    ylim(yl);
end
plot([0 0],[yl(1) yl(2)],'k--')
% gaps = findgaps(sig_ind);
% if ~isempty(gaps)
%     for g = 1:size(gaps,1)
%         gp = gaps(g,:);
%         gp(gp == 0) = [];
%         if length(gp) > 40
%             h = fill([min(gp) max(gp) max(gp) min(gp) min(gp)]-twin1,...
%                 [yl(1) yl(1) yl(2) yl(2) yl(1)],'k');
%             uistack(h,'down')
%             set(h,'facealpha',.25,'EdgeColor','None')
%         end
%     end
% end
xlim([-twin1 twin2]);
hold off
xlabel('Time from Fixation Start (ms)')
ylabel('Firing Rate (Hz)')
legend('out->in','out->out','Location','NorthWest')
set(gca,'Xtick',[-twin1 -twin1/2 0 twin1/2 twin1 3/2*twin1 twin2])
title(['Peak Firing Rate of ' num2str(pks,3) ' Hz at ' num2str(locs-twin1) ' ms'])
axis square

%% Plot Firing  Rate Map, Field Maps, and Raw Data

unit = cellfun(@(s) strfind(s,unit_name),unit_stats(1,:),'UniformOutput',false);
unit = find(cellfun(@(s) ~isempty(s), unit) == 1);

%---plot spike locations color coded by 1st half vs second half of session/by time---%
figure
make_spike_jittered_colored_plot(eyepos{unit},spike_times{unit},[1 1],1)
set(gca,'Xcolor','w')
set(gca,'Ycolor','w')


%%
binsize = 12; %pixels per bin spatial bin in either dimension 1/2 dva
filter_width = 6; %std of 2D guassian filter ~ 3 dva, could use 2 dva (4) as well
H = define_spatial_filter(filter_width);

imageX = 800;
imageY = 600;
%---Plot Rate Map for All Images---%
ratemap = get_firing_rate_map({eyepos{unit},spike_times{unit}},imageX,imageY,binsize,H,Fs,'all');

figure
h = imagesc(ratemap);
%                 hold on
%                 gry = 0.5*ones(size(ratemap,1),size(ratemap,2),3);
%                 h = imshow(gry);
%                 hold off
%                 set(h,'alphadata',isnan(ratemap));
set(h,'alphadata',~isnan(ratemap));
axis off
axis equal

clims = caxis;
max_fr = prctile(ratemap(:),97.5); % the ~97.5%-tile
caxis([clims(1) max_fr])
colormap('jet')
colorbar