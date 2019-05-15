function trl = trialfunclrchng(cfg)

% read the header
hdr = ft_read_header(cfg.dataset);
% read the events
event = ft_read_event(cfg.dataset);

eventt = event;
trial_start = 150;
trial_end = 151;
trial = {};
temp_event_code = [];
temp_event_time= [];
trial_count = 1;
while ~isempty(event);
    if event(1).value > 30000 %somehow bitshifted by 2^15
        if -1*(event(1).value-2^15+1) ~= trial_end;%so trial hasn't ended
            temp_event_code = [temp_event_code -1*(event(1).value-2^15+1)];
            temp_event_time = [temp_event_time event(1).sample];
            event(1) = [];
        else %if trial hasn't ended group all codes and times together into cell array
            temp_event_code = [temp_event_code -1*(event(1).value-2^15+1)];
            temp_event_time = [temp_event_time event(1).sample];
            event(1) = [];
            trial{1,trial_count} = temp_event_code;
            trial{2,trial_count} = temp_event_time;
            temp_event_time = [];
            temp_event_code = [];
            trial_count = trial_count+1;
        end
    else
        if event(1).value ~= trial_end;%so trial hasn't ended
            temp_event_code = [temp_event_code event(1).value];
            temp_event_time = [temp_event_time event(1).sample];
            event(1) = [];
        else %if trial hasn't ended group all codes and times together into cell array
            temp_event_code = [temp_event_code event(1).value];
            temp_event_time = [temp_event_time event(1).sample];
            event(1) = [];
            trial{1,trial_count} = temp_event_code;
            trial{2,trial_count} = temp_event_time;
            temp_event_time = [];
            temp_event_code = [];
            trial_count = trial_count+1;
        end
    end
end

event = eventt;

% careful when trying to find indeces and 5th event appears to have been
% coded wrong in some case it supposed to be 5000 + trial # by may just be
% trial # except for the 1st trial which may not have a pretrial in which
% the trial # is in the 3rd location
numrpt = size(trial,2);
valrptcnt = 0;
clear trl clrchgind
for rptlop = 1:numrpt
    if length(find(trial{1,rptlop} == 200)) ~=0 && length(find(trial{1,rptlop} == 151)) ~= 0 
        trlbegind = find(trial{1,rptlop} == 23); %only want the eye data when the dot is on
        trlendind = find(trial{1,rptlop} == 24); 
        if length(trlbegind) > 1
            trlbegind = trlbegind(2);
        end
        if length(trlendind) > 1
            trlendind = trlendind(2);
        end
        cndnumind = find(trial{1,rptlop} >= 1000 & trial{1,rptlop} <=2000);
        begtimdum = trial{2,rptlop}(trlbegind)-100;
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
