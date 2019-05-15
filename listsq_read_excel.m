function listsq_read_excel(data_dir,excel_file)
%function reads in session data stored in excel_file form and then stores
%the appropriate data into a .mat file;

% OUTPUT data file contains the following information
%
% 1) session_data: data organized by session in each cell
%   a) session_data{sess}.calibration_file: calibration file name
%   b) session_data{sess}.task1_file/task2_file: file name for first and
%   2nd task of the day (if there was a second)
%   c) session_data{sess}.task#: task save or item/cnd files
%   d) session_data{sess}.task#_unit_names: unit sames e.g. sig001a
%   e) session_data{sess}.task#_unit_waveforCount: number of waveforms
%   f) session_data{sess}.task#_unit_confidence: confidence that what is sorted is a real spike
%       i) ranges from 0-100% ...kind of arbitrary
%   g) session_data{sess}.task#_unit_cutQuality: quality of cluster cutting
%        ...usually based on seperability from noise cluster in PCA
%        space, but occassionally non-linear energy or peak-valley + PCA is
%        a better measure
%       i) 1: no seperability from noise
%       ii) 2: unit cluster is embedded in the noise but can see seperate cluster forming
%       iii) 3: unit cluster is juxtaposed to noise cluster
%       iv) 4: unit cluster is clearly seperated from noise cluster
%       v) 5: unit cluster is no where near noise cluster
%   h) session_data{sess}.task#_unit_MultiunitSeperability: quality of multiunit seperability
%       i) 0: data cleary show multiple units
%       ii) 1: >90% chance theres more than 1 unit
%       iii) 2: > 75% chance there is more than 1 unit
%       iv) 3: 50/50 whether it's 1 or multiple units
%       v) 4: >  90% chance that I think it is a single unit
%       vi)5: > 99% chance that I think it's a single unit based on the data
%   i) session_data{sess}.task#_LFP: 0 LFP is bad, 0 maybe keep LFP, 1 LFP is good for analysis

[~,~,raw] = xlsread(excel_file);


%find where each session starts using #### as the key
session_start = [];
for row = 1:size(raw,1)
    if strcmpi(raw{row,1},'####') || strcmpi(raw{row,1},'C####')
        session_start = [session_start, row];
    end
end

%remove empty sessions
if size(raw,1) == session_start(end);
    session_start(end) = [];
end

session_data = cell(1,length(session_start)); %storing all session data

for sess = 1:length(session_start);
    num_tasks = 0;%number of tasks per recording session usually 1 or 2
    
    start_row = session_start(sess)+1;
    if sess == length(session_start)
        end_row = size(raw,1);
    else
        end_row = session_start(sess+1);
    end
    
    sess_raw = raw(start_row:end_row,:);  %just the data for this session
    for row = 1:size(sess_raw,1);
        if ~isempty(strfind(lower(sess_raw{row,1}),'calibration'))
            session_data{sess}.calibration_file = string_trim(sess_raw{row,2});
        elseif strcmpi(sess_raw{row,1},'Filename')
            num_tasks = num_tasks+1;
            if num_tasks ==  1
                session_data{sess}.task1_file = string_trim(sess_raw{row+1,1});
                session_data{sess}.task1 = string_trim(sess_raw{row+1,6});
                session_data{sess}.location = [sess_raw{row+1,2} sess_raw{row+1,3}];
                session_data{sess}.subregion = sess_raw{row+1,4};
                session_data{sess}.task1_unit_names = []; %unit sames e.g. sig001a
                session_data{sess}.task1_unit_waveforCount = []; %number of waveforms
                session_data{sess}.task1_unit_confidence = []; %confidence that what is sorted is a real spike
                session_data{sess}.task1_unit_cutQuality = [];  %quality of cluster cutting using in PCA space from noise
                session_data{sess}.task1_unit_MultiunitSeperability = [];%quality of multiunit seperability
                session_data{sess}.task1_unit_comments = {}; %extra comments by person who sorted
                
            elseif num_tasks ==  2
                session_data{sess}.task2_file = string_trim(sess_raw{row+1,1});
                session_data{sess}.task2 = string_trim(sess_raw{row+1,6});
                session_data{sess}.task2 = string_trim(sess_raw{row+1,6});
                
                session_data{sess}.task2_unit_names = []; %unit sames e.g. sig001a
                session_data{sess}.task2_unit_waveforCount = []; %number of waveforms
                session_data{sess}.task2_unit_confidence = []; %confidence that what is sorted is a real spike
                session_data{sess}.task2_unit_cutQuality = [];  %quality of cluster cutting using in PCA space from noise
                session_data{sess}.task2_unit_MultiunitSeperability = [];%quality of multiunit seperability
                session_data{sess}.task2_unit_comments = {}; %extra comments by person who sorted
            else
                disp('More than 2 tasks run')
            end
        elseif ~isempty(strfind(upper(sess_raw{row,1}),'LFP'))
            if num_tasks ==  1
                session_data{sess}.task1_LFP = cell2mat(sess_raw(row+1,2:5));
            elseif num_tasks ==  2
                session_data{sess}.task2_LFP = cell2mat(sess_raw(row+1,2:5));
            end
        elseif ~isempty(strfind(sess_raw{row,1},'sig'))
            if num_tasks ==  1
                session_data{sess}.task1_unit_names = [session_data{sess}.task1_unit_names sess_raw(row,1)];
                session_data{sess}.task1_unit_waveforCount= [session_data{sess}.task1_unit_waveforCount sess_raw{row,2}];
                if sess_raw{row,2} == 0 %no waveforms
                    session_data{sess}.task1_unit_confidence = [session_data{sess}.task1_unit_confidence NaN];
                    session_data{sess}.task1_unit_cutQuality = [session_data{sess}.task1_unit_cutQuality NaN];
                    session_data{sess}.task1_unit_MultiunitSeperability = [session_data{sess}.task1_unit_MultiunitSeperability NaN];
                    session_data{sess}.task1_unit_comments = [session_data{sess}.task1_unit_comments  {NaN}];
                else
                    session_data{sess}.task1_unit_confidence = [session_data{sess}.task1_unit_confidence sess_raw{row,3}];
                    session_data{sess}.task1_unit_cutQuality = [session_data{sess}.task1_unit_cutQuality sess_raw{row,4}];
                    session_data{sess}.task1_unit_MultiunitSeperability = [session_data{sess}.task1_unit_MultiunitSeperability sess_raw{row,5}];
                    session_data{sess}.task1_unit_comments = [session_data{sess}.task1_unit_comments  sess_raw(row,6)];
                end
            elseif num_tasks ==  2
                session_data{sess}.task2_unit_names = [session_data{sess}.task2_unit_names sess_raw(row,1)];
                session_data{sess}.task2_unit_waveforCount= [session_data{sess}.task2_unit_waveforCount sess_raw{row,2}];
                if sess_raw{row,2} == 0 %no waveforms
                    session_data{sess}.task2_unit_confidence = [session_data{sess}.task2_unit_confidence NaN];
                    session_data{sess}.task2_unit_cutQuality = [session_data{sess}.task2_unit_cutQuality NaN];
                    session_data{sess}.task2_unit_MultiunitSeperability = [session_data{sess}.task2_unit_MultiunitSeperability NaN];
                    session_data{sess}.task2_unit_comments = [session_data{sess}.task2_unit_comments  {NaN}];
                else
                    session_data{sess}.task2_unit_confidence = [session_data{sess}.task2_unit_confidence sess_raw{row,3}];
                    session_data{sess}.task2_unit_cutQuality = [session_data{sess}.task2_unit_cutQuality sess_raw{row,4}];
                    session_data{sess}.task2_unit_MultiunitSeperability = [session_data{sess}.task2_unit_MultiunitSeperability sess_raw{row,5}];
                    session_data{sess}.task2_unit_comments = [session_data{sess}.task2_unit_comments  sess_raw(row,6)];
                end
            end
        end
    end
end

slashes = strfind(excel_file,'\');
monkey = excel_file(slashes(end-1)+1:slashes(end)-1);
save([data_dir 'Across_Session_Unit_Data_' monkey '.mat'],'session_data')
end

function output = string_trim(input)
if ~isnan(input)
    output = strtrim(input); %remove white space
else
    output = input;
end
end
