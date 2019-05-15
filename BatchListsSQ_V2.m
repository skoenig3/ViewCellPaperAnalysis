% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task.
% code runs all the other code preprocess then process recording data
clar

task = 'ListSQ';
set(0,'DefaultFigureVisible','OFF');
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = 'P:\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PW Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\PW Figures\';
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = 'P:\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ViewCellPaperAnalyses\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
        %listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import and Pre-Process Recording Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         ImportListSQRecordingDataV2(data_dir,figure_dir,session_data{session})
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Correct Event Codes/Spike-LFP/Eye Data Misalignment---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session =1:length(session_data)
    %         disp(['Running Lag Check on session#' num2str(session) ' monkey#' num2str(monkey)])
    %         SpikeLFPLagCorrection(data_dir,figure_dir,session_data{session})
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Plot Waveforms and Rasters to Determine Firing Stability---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     addpath('C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PlexonSDK_v1_6\')
    %             for session = 1:length(session_data)
    %                 disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %                 make_rasters_and_plot_waveformsV2(data_dir,figure_dir,session_data{session},task)
    %             end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Spatially Modulated---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session = 1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         spatial_analysisV2(data_dir,figure_dir,session_data{session},task)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%---Determine Place Cell Reliablity Across Fixations in Field vs out Field---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for session =1:length(session_data)
    %         disp(['#' num2str(session) ' ' session_data{session}.task1_file])
    %         Place_Cell_Fixation_Analysis(data_dir,figure_dir,session_data{session},predict_rt)
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Determine if Neurons are Fixation/Saccade Modulated in List Task---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1:length(session_data)
        List_Fixation_Analysis(data_dir,figure_dir,session_data{session})
        List_Saccade_AnalysisV2(data_dir,figure_dir,session_data{session}) %not rechecked for bugs
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Modulated by Saccade Direction and/or Amplitude---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1:length(session_data)
        disp(['#' num2str(session) ' ' session_data{session}.task1_file])
        ListSQ_Saccade_Direction_and_Amplitude_Analysis_Future_and_Past(data_dir,figure_dir,session_data{session})
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Neurons are Visually reponsive/like Cross hair---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1:length(session_data)
        Visual_Response_AnalysisV2(data_dir,figure_dir,session_data{session});
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Determine if Visually Responsive Neurons Show Novel/Repeat Differences---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for session = 1 :length(session_data)
        Visual_Response_Memory(data_dir,figure_dir,session_data{session});
    end
    
end
set(0,'DefaultFigureVisible','ON');