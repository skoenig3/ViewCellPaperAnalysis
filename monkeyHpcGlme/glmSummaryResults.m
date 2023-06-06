% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task.
% code runs all the other code preprocess then process recording data
clar

task = 'ListSQ';
monkeys = {'Vivian','Tobii'};
figure_dir = 'D:\MATLAB\ViewCellPaperAnalysis\glmAnalayis\';
for monk = 2:-1:1
    monkey = monkeys{monk};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if strcmpi(monkey,'Vivian')
        data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\PW Recording Files\';
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        

    elseif strcmpi(monkey,'Tobii')
        data_dir = 'D:\MATLAB\ViewCellPaperAnalysis\TO Recording Files\';
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine if Neurons are Fixation/Saccade Modulated in List Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1:length(session_data)
         getDataForGlmeListSQ(monk,session,data_dir,session_data{session})
    end
    

end