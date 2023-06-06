%Script creates and runs monkey HPC eye-movement firingrate model using  a
%Linear/Non-linear GLM Poisson Model based on rodent code Hardcastle et al.
%2017. Some insparation from Michael's paper too (You et al., 2020).
%A few major differences:
%   1) Using Monte Carlo Cross Validation (MCCV) instead of k-fold because
%   less data and data structure is biased over time by block and such
%   2) Poisson data is smoothed and represent firing rate within fixation &
%   previous saccade as eye movement seem like a natural binning process.
%Written by Seth Koenig 7/24/2022

clar


%---Processing Parameters---%
mainDataDir = 'E:\monkeyHpcData\';
figureDir = 'E:\monkeyHpcData\Figures\';
monkeyNames = {'PW','TO'};

%monte carlo cross-validation parameters
mCCVParms = [];
mCCVParms.shuffs = 100;
mCCVParms.dataSplit = 0.70;%70/30 train test
mCCVParms.variablesNames = {'P','D','E'}; %position, saccade direction, and eye time/phase
mCCVParms.typeParams = {'2d','1dcirc','1dcirc'}; %needed for roughness regularization
mCCVParms.regularizationType = 'roughness'; %supporsts L1, L2, and roughness (preferred)
mCCVParms.regularizationWeights = [10 9 7]; %for roughness regularization, somewhat arbitrary @ ~2*sqrt(2)

%bin size parameters
binSizeParams = [];
binSizeParams.dimX = 800; %in pixels
binSizeParams.dimY = 600; %in pixels
binSizeParams.numSpatialBins = 25;%~1 dva, square so x/y not equal
binSizeParams.numDirBins = 20;%18 degrees per bin or 18 bins
binSizeParams.numTimesBins = 12;%~480 ms or so


%%
%---Run the Model Data---%
allBestModels = [];
allMeanParams = [];
allUnitNames = [];
allMonkeyFixationData = [];
allMonkeyFixationLockedFiringBinned = [];
allMonkeyFixationLockedFiring = [];
allMonkeyFixationLockedFiringRaw = [];
allSpatialCorr = [];
allSigSpatialScores = [];
for mN = 2:-1:1
    %get list of files/session data
    thisMonkeyDir = [mainDataDir monkeyNames{mN} '_glmeData2' filesep];
    fileList = dir(fullfile(thisMonkeyDir,'*.mat'));
    
    %load each session and then run the models
    for fL = 1:length(fileList)
        load([thisMonkeyDir fileList(fL).name])
        
        %store unit names
        for unit = 1:length(allFixationLockedFiringBinned)
            if ~isempty(allFixationLockedFiringBinned{unit})
                allUnitNames = [allUnitNames {[fileList(fL).name(1:8)  unit_stats{1,unit}]}];
            else
                allUnitNames = [allUnitNames {' '}];
            end
        end
        
        %then go through each unit
        [bestModels,meanParams] = runLNGlmeMonkeyHPC(figureDir,fileList(fL).name(1:8),...
            allFiringRateMaps,allFixationLockedFiringBinned,mCCVParms,binSizeParams);
        allBestModels = [allBestModels bestModels];
        allMeanParams = [allMeanParams meanParams];
        
        %store other data
        allMonkeyFixationData = [allMonkeyFixationData allFixationData];
        allMonkeyFixationLockedFiringBinned = [allMonkeyFixationLockedFiringBinned allFixationLockedFiringBinned];
        allMonkeyFixationLockedFiring = [allMonkeyFixationLockedFiring allFixationLockedFiring];
        allMonkeyFixationLockedFiringRaw = [allMonkeyFixationLockedFiringRaw allFixationLockedFiringRaw];
        allSpatialCorr = [allSpatialCorr allSpatialCorr];
        allSigSpatialScores = [allSigSpatialScores sigSpatialScores];
    end
end


%---Save the Processed Data---%
save([mainDataDir 'LnSummaryData.mat'],...
    'allBestModels','allBestModels','allMeanParams','allUnitNames',...
    'allMonkeyFixationData','allMonkeyFixationLockedFiringBinned',...
    'allMonkeyFixationLockedFiring','allMonkeyFixationLockedFiringRaw',...
    'allSpatialCorr','allSigSpatialScores','mCCVParms','binSizeParams')

%%
%---Get Summary Unit Counts---%
%get model names and parameters
[A,modelTypes,modelNames,numModels,numModelParameters] = ...
    buildModelDataMonkeyHpc(cell(1,length(mCCVParms.variablesNames)),mCCVParms.variablesNames);
modelTypes2 = cell2mat(modelTypes);

%get overal counts
numUnits = sum(~isnan(allBestModels));%with enough spikes
numSigUnits = sum(allBestModels > 0);
numNotEnoughSpikes = sum(allBestModels == -99);
numNonSigUnits = sum(allBestModels == -1);

%get sig counts by model type
sigCounts = NaN(1,length(modelNames));
for model = 1:length(modelNames)
    sigCounts(model) = sum(allBestModels == model);
end

%num units with each type
numLnSpatialUnits = sum(sigCounts(1,modelTypes2(:,1) == 1)); %any spatial
numDirectionUnits = sum(sigCounts(1,modelTypes2(:,2) == 1)); %any direction
numEyeUnits = sum(sigCounts(1,modelTypes2(:,3) == 1));%any eye

%determine congruency with "tuning curve" anlaysis
spatialModelsNumbers = find(modelTypes2(:,1) == 1);
nonSpatialDirectionModelsNumbers = find(modelTypes2(:,1) == 0 & modelTypes2(:,2) == 1);
sumSigTuningCurve = sum(allSigSpatialScores == 1);

%wasnt spatial before but now is
newlySigModel = allBestModels(1,allBestModels > 0 & allSigSpatialScores == 0);
newlySpatialModel = 0;
for sM = 1:length(spatialModelsNumbers)
    newlySpatialModel = newlySpatialModel + sum(newlySigModel == spatialModelsNumbers(sM));
end

%still spatial
hasAnySigModel = allBestModels(1,allBestModels > 0 & allSigSpatialScores == 1); %and TC is sig spatial 
stillSigSpatial = 0;
for sM = 1:length(spatialModelsNumbers)
    stillSigSpatial = stillSigSpatial + sum(hasAnySigModel == spatialModelsNumbers(sM));
end

%no longer spatial but is directional
hasAnySigModel = allBestModels(1,allBestModels > 0 & allSigSpatialScores == 1); %and TC is sig spatial 
nowDirectional = 0;
for nSDMN = 1:length(nonSpatialDirectionModelsNumbers)
    nowDirectional = nowDirectional + sum(hasAnySigModel == nonSpatialDirectionModelsNumbers(nSDMN));
end

%LN eye only
hasAnySigModel = allBestModels(1,allBestModels > 0 & allSigSpatialScores == 1); %and TC is sig spatial 
eyeOnly = sum(hasAnySigModel == 7);

%other categories
couldntModel = sum(allSigSpatialScores == 1 & allBestModels < 0);
noLongerSigSpatial = sum(allSigSpatialScores == 1 & isnan(allBestModels));

if stillSigSpatial + couldntModel + noLongerSigSpatial + nowDirectional  + eyeOnly ~= sumSigTuningCurve
   error('counts are off') 
end


%%
%---Plot Summary Count Data---%
figure
subplot(2,2,1)
bar(sigCounts')
xticks(1:length(modelNames))
xticklabels(modelNames)
ylabel('Number of Modeled Units')
title('Distribution of Signiciant Units')

subplot(2,2,2)
bar([numSigUnits numNonSigUnits numNotEnoughSpikes]/numUnits)
text(1,numSigUnits/numUnits+0.05,['n = ' num2str(numSigUnits)])
text(2,numNonSigUnits/numUnits+0.05,['n = ' num2str(numNonSigUnits)])
text(3,numNotEnoughSpikes/numUnits+0.05,['n = ' num2str(numNotEnoughSpikes)])
xticks(1:3)
ylim([0 0.6])
ylabel('Prop All Units')
xticklabels({'Sig Units','Non-Sig Units','Too Few Spikes'})
title('Distribution of Units Modeled')

subplot(2,2,3)
bar([numLnSpatialUnits numDirectionUnits numEyeUnits]/numUnits)
hold on
text(1,numLnSpatialUnits/numUnits+0.05,['n = ' num2str(numLnSpatialUnits)])
text(2,numDirectionUnits/numUnits+0.05,['n = ' num2str(numDirectionUnits)])
text(3,numEyeUnits/numUnits+0.05,['n = ' num2str(numEyeUnits)])
hold off
xticks(1:3)
ylim([0 0.6])
xticklabels({'Any Spatial','Any Direction','Any Eye'})
title('LN Model: Prop Units that Have Any Encoding Type')

subplot(2,2,4)
bar([sumSigTuningCurve numLnSpatialUnits stillSigSpatial newlySpatialModel noLongerSigSpatial couldntModel nowDirectional eyeOnly])
xticks(1:8)
xticklabels({'Sig TC','Sig LN','Sig TC & LN','Sig LN Only','Sig TC Only','TC but FR too Low for LN','LN Dir Not Space','LN Eye Only'})
ylabel('Num Units')
title('Breakdown of Spatial Units Across Anlayses')

sgtitle('Monkey HPC Linear Non-linear GLME Model Summary Results')

%save figure as .png and .fig
save_and_close_fig(mainDataDir,'summaryLNModelResults')