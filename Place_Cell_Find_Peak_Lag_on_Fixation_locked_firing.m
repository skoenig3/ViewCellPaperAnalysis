data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';

w = warning ('off','all');

a = what(data_dir);
a = a.mat;

all_lags = [];
min_all_lags = []; 
n_cells = 0;
n_seq_in_out = 0; 
p_seq_in_out = 0;
for file = 1:length(a)
    if ~isempty(strfind(a{file},'-View_Cell_Analysis.mat'))
         load([data_dir a{file}],'list_fixation_locked_firing','in_out','twin',...
             'smval','multiunit','task_file','seq_p_in_out','sequences_inside')
    else
        continue
    end
    
    for unit = 1:length(list_fixation_locked_firing)
        if multiunit(unit) 
            continue
        end
        if ~isempty(list_fixation_locked_firing{unit})
            n_cells = n_cells+1;
            if ~all(isnan(seq_p_in_out(unit,:)))
                n_seq_in_out = n_seq_in_out+1;
                pvals = seq_p_in_out(unit,50:end);%after fixation onset
                if sum(pvals < 0.05) > 2
                    p_seq_in_out = p_seq_in_out+1;
                end
            end
            
            t = -twin:twin-1;
            figure
            [~,~,~,y] = dofill(t,list_fixation_locked_firing{unit}(in_out{unit} == 1,:),'blue',1,smval);%inside
            [peakLoc] = peakfinder(y);
            peaks = y(peakLoc);
            hold on
            dofill(t,list_fixation_locked_firing{unit}(in_out{unit} == 0,:),'red',1,smval);%outside
            for p = 1:length(peaks)
                plot(peakLoc-twin,peaks,'*k')
            end
            hold off
            xlabel('Time from Fixation Onset (ms)')
            lag = peakLoc(peaks == max(peaks))-twin;
            lag = lag(1);
            peakLoc = peakLoc-twin;
            peakLoc(peakLoc < 0) = [];
            if ~isempty(peakLoc)
                min_lag = peakLoc(1);
                if min_lag < 0
                   disp('what?')
                end
            else
                min_lag = NaN;
            end
            ylabel('Firing Rate (Hz)')
            legend('Fix Inside','Fix Outside','Location','NorthEastOutside')
            title(sprintf(['n_{in} = ' num2str(sum(in_out{unit} == 1)) ', n_{out} = ' num2str(sum(in_out{unit} == 0)) ...
                ' Peak Lag = ' num2str(lag) ' Min Lag = ' num2str(min_lag)]))
            xlim([-twin twin]);
            
            all_lags = [all_lags lag];
            min_all_lags = [min_all_lags min_lag]; 
            
%             pause(3)
            close

        end
    end
    
end