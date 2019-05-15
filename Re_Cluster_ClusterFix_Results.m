%written by Seth Konig August 2016
%Cluster Fix seems to have done a good job dividing scan paths into fixations
%and saccades, but seems to have missed a few of the smaller saccades that divide
%long "fixation periods" into multiple fixations. We want to see if we can
%detect these small saccades but re-analyizing the Cluster Fix results for
%longer fixations. Just want to do the best that we can especially since
%certain analyses (e.g. fixation analysis for visual place cells) depends
%on good Cluster Fix results. We want that last few %.
%
%There seemed to be very few issues ~0.3% in fixations <
%500 ms in duration so going to look at the ones longer than 400 ms since
%want a buffer on that. Also, there are more shorter fixations than longer
%fixations during free viewing so any very small % will probably not affect
%the results as much. Estimated missed 3.87% of longer fixations/saccades or
% 1.1% of all fixations/saccades.
%
%Also, going to ignore any potential saccades that are less than 2 dva
%since these are really hard to detect and by visual inspection not as
%reliable similar to what we found with the Cluster Fix paper. These
%mini-saccades at least should have less of an affect on the results since
%eyes don't really move. Small saccades are ignored for a lot of analyses
%too.

max_dur = 400; %maximum duration we will not re-analyze fixations for
twin = 25; %buffer duration around "fixation" of interest

dispersion_threshold = 48;%48 pixels = 2 dva. If "fixation" is not very tight ...
%aka covers more than 2 pixels in space than look at it more carefully
%because it could really be more than 1 fixation

min_saccade_amp = 48; %48 pixels = 2 dva. Minimum saccade amplitude to
%count as a saccade kind of arbitrary but based on what I see from the data.

total_fixations = 0; %get count so can get accuracy
total_missed = 0;%number of newly detected fixations
total_long_fixations = 0; %number of fixaitons longer than max_dur;

task = 'ListSQ';
tic
for monkey = 1%:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resored Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 156;%156 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 138;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working on importing the data
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import and Pre-Process Recording Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for session = 1%:length(session_data)
        %get important task related data
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,sorting_quality,~]=...
            get_task_data(session_data{session},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'fixationstats');
        
        
        
        for t = 1:length(fixationstats)
            xy = fixationstats{t}.XY; %eye position
            fixations = fixationstats{t}.fixations; %mean fixation location
            fixationtimes = fixationstats{t}.fixationtimes; %time at which fixations occur
            
            if isempty(fixationtimes)
                continue
            end
            
            total_fixations = total_fixations+size(fixations,2);
            
            %variables to store newly found fixation and saccade times
            new_saccadetimes = [];
            new_fixationtimes = [];
            fix_index = [];
            
            %find fixations that are longer than the maximum duration
            dur = fixationtimes(2,:)-fixationtimes(1,:)+1; %fixation duration
            possibly_too_long = find(dur >= max_dur);
            
            total_long_fixations = total_long_fixations+length(possibly_too_long);
            
            
            for ptl = 1:length(possibly_too_long)
                fixt = fixationtimes(1,possibly_too_long(ptl)):fixationtimes(2,possibly_too_long(ptl));
                
                %can't really look at fixations on the edges and hopefully they
                %don't matter as much since in ITI or reward
                if fixt(1) < twin || fixt(end) > length(xy)-twin
                    continue
                end
                
                x_pos = xy(1,fixt);
                y_pos = xy(2,fixt);
                
                %Calculate how much the "fixation" spans to determine if it needs
                %to be re-analyzed. A single fixation should have tight spread.
                x_disp = max(x_pos)-min(x_pos);
                y_disp = max(y_pos)-min(y_pos);
                
                %if dispersion is greater than threshold re-analyze for more fixations
                if (x_disp > dispersion_threshold) || (y_disp > dispersion_threshold)
                    %                     figure
                    %                     plot(x_pos,y_pos,'k')
                    
                    %add saccade buffer so can properly run Cluster Fix
                    x_pos = xy(1,fixt(1)-25:fixt(end)+twin);
                    y_pos = xy(2,fixt(1)-25:fixt(end)+twin);
                    
                    %Re-run ClusterFix on fixation
                    eyedat = [x_pos;y_pos];
                    [fixationstats1] = ClusterFix_Plexon(eyedat);
                    
                    if size(fixationstats1.fixations,2) > 1 %means detected at least 2 saccades where there was 1 before
                        XY =  fixationstats1.XY(:,twin+1:end-twin);
                        
                        
                        %remove times that were part of saccade buffer
                        sactimes = fixationstats1.saccadetimes-twin;
                        sactimes(:,sactimes(1,:) < 1) = []; %saccade started before buffer started
                        sactimes(:,sactimes(1,:) > length(XY)) = []; %saccade started after buffer started
                        sactimes(:,sactimes(2,:) > length(XY)) = []; %saccade end ended after buffer so set to buffer
                        
                        
                        %find saccades that are really too small to be detected by
                        %our eye tracking system and therefore could be noise
                        too_small = zeros(1,size(sactimes,2));
                        for s = 1:size(sactimes,2)
                            x_amp = XY(1,sactimes(2,s))-XY(1,sactimes(1,s));
                            y_amp = XY(2,sactimes(2,s))-XY(2,sactimes(1,s));
                            amp = sqrt(x_amp^2+y_amp^2);
                            if amp < min_saccade_amp
                                too_small(s) = 1;
                            end
                            
                        end
                        
                        %Correct fixation times for saccade buffer will ignore any
                        %changes to start and end b/c I'm not sure what the
                        %appropriate window is for this potentially grouped fixation
                        %because we are not gauranteed that the neighboring saccades
                        %are good examples
                        fixtimes = fixationstats1.fixationtimes-twin;
                        if fixtimes(1,1) < 1
                            %it possible that fixation would start in this window before it started
                            %in larger window in Cluster Fix but going to ignore
                            %                     disp('Fixation Started earlier than originally thought...')
                            fixtimes(1,1) = 1;
                        elseif fixtimes(1,1) ~= 1;
                            %it possible that fixation would start in this window later than it started
                            %in larger window in Cluster Fix but going to ignore
                            %                     disp('Fixation Started later than originally thought...')
                            fixtimes(1,1) = 1;
                        end
                        if fixtimes(2,end) > length(XY)
                            %it possible that fixation would end in this window later than it ended
                            %in larger window in Cluster Fix but going to ignore
                            
                            %                     disp('Fixation ended later than originally thought...')
                            fixtimes(2,end) = length(XY);
                        elseif fixtimes(2,end) ~= length(XY)
                            %it possible that fixation would end in this window before it ended
                            %in larger window in Cluster Fix but going to ignore
                            %                     disp('Fixation ended earlier than originally thought...')
                            fixtimes(2,end) = length(XY);
                        end
                        
                        
                        %put all fixation time points into vector for fixations
                        fixationtimes2 = zeros(1,length(XY));
                        for f = 1:size(fixtimes,2);
                            fixationtimes2(fixtimes(1,f):fixtimes(2,f)) = 1;
                        end
                        
                        %now add too small "saccades" to fixation time vector to
                        %concantenate with neighboring fixations
                        if any(too_small);
                            small = find(too_small);
                            for s = 1:length(small);
                                fixationtimes2(sactimes(1,small(s)):sactimes(2,small(s))) = 1;
                            end
                            
                            %return to standard format
                            ft2 = find(fixationtimes2);
                            broken_ind=findgaps(ft2);
                            fixtimes = NaN(2,size(broken_ind,1));
                            for f = 1:size(broken_ind,1);
                                bi = broken_ind(f,:);
                                bi(bi == 0) = [];
                                fixtimes(:,f) = [bi(1); bi(end)];
                            end
                        end
                        
                        sac_dur = sactimes(2,:)-sactimes(1,:);
                        sactimes(:,sac_dur <10) = [];
                        
                        
                        if size(fixtimes,2) > 1 %have newly detected saccades/fixations
                            
                            %plot newly discovered fixations in red on original eye trace in black
                            total_missed = total_missed+size(fixtimes,2)-1;
                            %                             hold on
                            %                             for f = 1:size(fixtimes,2)
                            %                                 plot(XY(1,fixtimes(1,f):fixtimes(2,f)),XY(2,fixtimes(1,f):fixtimes(2,f)),'r');
                            %                             end
                            
                            %Current times relative to this fixation so get times
                            %relative to task start using start time of original
                            %fixation
                            fixtimes = fixtimes+ fixationtimes(1,possibly_too_long(ptl))-1; %-1 since 1st index is the same time/overlaps
                            sactimes = NaN(2,size(fixtimes,2)-1);
                            for f = 1:size(fixtimes,2)-1;
                                sactimes(:,f) = [fixtimes(2,f)+1;fixtimes(1,f+1)-1];
                            end
                            
                            %store variables across fixations within trial
                            new_fixationtimes = [new_fixationtimes {fixtimes}];
                            new_saccadetimes = [new_saccadetimes {sactimes}];
                            fix_index = [fix_index possibly_too_long(ptl)];
                        end
                    end
                end
            end
            
            %now put into original fixationstats into trial across the fixations in
            %this trial
            if ~isempty(new_saccadetimes)
                for f = length(fix_index):-1:1 %go in reverse order so easier
                    new_fixations = NaN(2,size(new_fixationtimes{f},2));
                    for ff = 1:size(new_fixationtimes{f},2)
                        new_fixations(:,ff) = mean(xy(:,new_fixationtimes{f}(1,ff):new_fixationtimes{f}(2,ff)),2);
                    end
                    
                    fixationstats{t}.fixations = [fixationstats{t}.fixations(:,1:fix_index(f)-1) ...
                        new_fixations fixationstats{t}.fixations(:,fix_index(f)+1:end)];
                    fixationstats{t}.fixationtimes = [fixationstats{t}.fixationtimes(:,1:fix_index(f)-1) ...
                        new_fixationtimes{f} fixationstats{t}.fixationtimes(:,fix_index(f)+1:end)];
                    
                    %just squeeze the new saccadetimes in
                    original_saccadetimes = fixationstats{t}.saccadetimes;
                    sac_order = find((new_saccadetimes{f}(1) > original_saccadetimes(1,:)));
                    if isempty(sac_order)
                        fixationstats{t}.saccadetimes = [new_saccadetimes{f} fixationstats{t}.saccadetimes];
                    else
                        fixationstats{t}.saccadetimes = [fixationstats{t}.saccadetimes(:,1:sac_order(end)) ...
                            new_saccadetimes{f} fixationstats{t}.saccadetimes(:,sac_order(end)+1:end)];
                    end
                end
                
                
                xy = fixationstats{t}.XY; %eye position
                fixations = fixationstats{t}.fixations; %mean fixation location
                fixationtimes = fixationstats{t}.fixationtimes; %time at which fixations occur
                saccadetimes = fixationstats{t}.saccadetimes;
                
                %---Plot Re-Analyzed Cluster Fix Results for the whole scan path---%
                %                 figure
                %                 hold on
                %                 plot(xy(1,:),xy(2,:),'k')
                %                 for f = 1:size(fixationtimes,2)
                %                     plot(xy(1,fixationtimes(1,f):fixationtimes(2,f)),xy(2,fixationtimes(1,f):fixationtimes(2,f)),'g')
                %                     plot(fixations(1,f),fixations(2,f),'kd','markersize',8)
                %                 end
                %                 for f = 1:size(saccadetimes,2)
                %                     plot(xy(1,saccadetimes(1,f):saccadetimes(2,f)),xy(2,saccadetimes(1,f):saccadetimes(2,f)),'m')
                %                 end
                %                 xlim([0 800])
                %                 ylim([0 600])
                %
                %                 figure
                %                 hold on
                %                 plot(xy(1,:),'k')
                %                 plot(xy(2,:),'k')
                %                 for f = 1:size(fixationtimes,2)
                %                     plot(fixationtimes(1,f):fixationtimes(2,f),xy(1,fixationtimes(1,f):fixationtimes(2,f)),'g')
                %                     plot(fixationtimes(1,f):fixationtimes(2,f),xy(2,fixationtimes(1,f):fixationtimes(2,f)),'g')
                %                 end
                %                 for f = 1:size(saccadetimes,2)
                %                     plot(saccadetimes(1,f):saccadetimes(2,f),xy(1,saccadetimes(1,f):saccadetimes(2,f)),'m')
                %                     plot(saccadetimes(1,f):saccadetimes(2,f),xy(2,saccadetimes(1,f):saccadetimes(2,f)),'m')
                %                 end
                %                 ylim([0 800])
                %
                %                 close all
            end
        end
        save([data_dir task_file(1:end-11) '-preprocessed.mat'],'-append','fixationstats'); %resave the fixationstats to the file
    end
end
tm = toc;
emailme(['ReClusterFixingtook ' num2str(tm)])
emailme(['Total fixations missed: ' num2str(total_missed) ' total fixations: ' num2str(total_fixations)])