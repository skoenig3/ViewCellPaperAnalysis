function trl = trialfunListSQ(cfg)
% Cuts out all of the trial data from one session.
% Structure of a trial in event codes:
% 15 start pre-trial
% 16 end pre-trial
% 150 TRIAL_START
% 151 Trial end
% 100 start eye data
% 5000+x trial number
% 1000+x condition number
% 500+x block number


% read the header
hdr = ft_read_header(cfg.dataset);
% read the events
event = ft_read_event(cfg.dataset);

% going to add codes to the end of trials since there is ~100 ms break in
% between for cortex
start_pause = 19;
end_pause = 20;
event = struct2cell(event);
event = cell2mat(event(1:2,:));
if abs(event(2,1)) > 30000 %16 bit shift
    event(2,:) = -1*(event(2,:)-2^15+1);
end
trial_start = 15;
trial_end = 151;
trial = {};
trial_count = 1;
starts = find(event(2,:) == trial_start);
ends = find(event(2,:) == trial_end);
if starts(1) > ends(1); %must have started recording after the task started
    if event(2,1) == 150 %skipped pre-trial for whatever reason
        starts =[1 starts]; %append first event to starts
    else
        ends(1) = []; %remove ends since can't use trial
    end
end
if length(starts)-length(ends) == 1;
    starts(end) = [];
elseif  length(starts)-length(ends) > 1
   error('error something wrong with the encodes. More than 1 difference in when trial starts and trials end');
end
for t = 1:length(starts);
    trial{1,t} = event(2,starts(t):ends(t)); %event codes
    trial{2,t} = event(1,starts(t):ends(t)); %time of event
end

% careful when trying to find indeces and 5th event appears to have been
% codede wrong in some case it supposed to be 5000 + trial # by may just be
% trial # except for the 1st trial which may not have a pretrial in which
% the trial # is in the 3rd location
numrpt = size(trial,2);
valrptcnt = 0;
clear trl clrchgind
for rptlop = 1:numrpt
    if (trial{1,rptlop}(find(trial{1,rptlop}>500,1,'last'))-500) > 1 %1st block is color change for offset correction so ignore
        trlbegind = find(trial{1,rptlop} ==  15); % start at pretrial since I need this data for some things
        trlendind = find(trial{1,rptlop} == 151); % end at end trial
        if length( trlbegind) > 1
            trlbegind = trlbegind(1);
            trlendind = trlendind(1);
        end
        cndnumind = find(trial{1,rptlop} >= 1000 & trial{1,rptlop} <=2000);
        begtimdum = trial{2,rptlop}(trlbegind);
        endtimdum = trial{2,rptlop}(trlendind);
        if endtimdum > begtimdum
            valrptcnt = valrptcnt + 1;
            clrchgind(valrptcnt)=rptlop;
            trl(valrptcnt).begsmpind = begtimdum;
            trl(valrptcnt).endsmpind = endtimdum;
            trl(valrptcnt).cnd = trial{1,rptlop}(cndnumind);
            trl(valrptcnt).allval = trial{1,rptlop};
            trl(valrptcnt).alltim = trial{2,rptlop};
            trl(valrptcnt).event = rptlop;
        end
    end
end

for t = 1:length(trl)
    if t == length(trl)
        trl(t).endsmpind = trl(t).endsmpind+100;
        trl(t).allval = [trl(t).allval start_pause end_pause];
        trl(t).alltim = [trl(t).alltim trl(t).alltim(end)+1 trl(t).alltim(end)+100];
    else
        next_trialstart = trl(t+1).alltim(1);
        trl(t).endsmpind = next_trialstart-1;
        trl(t).allval = [trl(t).allval start_pause end_pause];
        trl(t).alltim = [trl(t).alltim trl(t).alltim(end)+1 next_trialstart-1];
    end
end
