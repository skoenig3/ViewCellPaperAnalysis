%Script analyzes population of linear-nonlinear models and builds decoders
%from them.
%Written by Seth Konig 06/01/2023


clear, clc, close all

%random seed
rng(06012023,'twister')


%---What and Where---%
mainDataDir = 'E:\monkeyHpcData\';
monkeyNames = {'PW','TO'};
savedLnModelResults = 'LnSummaryData.mat';
plotFlag = false;

%decoding over time window length, should be smval
erpParameters = getERPParameters(1000,400,200);

%shuffling count for decoding error since uneven unit counts
numShuffs = 1000;
tBin = 0.05;%50 ms time bins

%%
%---Load the Data---%
load([mainDataDir savedLnModelResults])


%%
%---Run the Model Data---%
allUnitNames2 = [];
allFiringRateMapsSmoothed = [];
allMonkeyFixationDataSmoothed = [];
allMonkeyFixationLockedFiringSmoothed = [];
for mN = 2:-1:1
    %get list of files/session data
    thisMonkeyDir = [mainDataDir monkeyNames{mN} '_glmeData2' filesep];
    fileList = dir(fullfile(thisMonkeyDir,'*.mat'));

    %load each session and then run the models
    for fL = 1:length(fileList)
        load([thisMonkeyDir fileList(fL).name])

        %store smoothed data
        for unit = 1:size(unit_stats,2)
            allUnitNames2 = [allUnitNames2; {[fileList(fL).name(1:8) unit_stats{1,unit}]}];
        end
        allFiringRateMapsSmoothed = [allFiringRateMapsSmoothed allFiringRateMaps];
        allMonkeyFixationDataSmoothed = [allMonkeyFixationDataSmoothed allFixationData];
        allMonkeyFixationLockedFiringSmoothed = [allMonkeyFixationLockedFiringSmoothed allFixationLockedFiring];
    end
end

%check length of data matches
if length(allUnitNames) ~= length(allUnitNames2)
    error('wtf')
end


%check same units
for unit = 1:length(allUnitNames2)
    if isempty(allUnitNames{unit}) || strcmpi(allUnitNames{unit},' ')
        %means did not have high enough firing rate to test
        continue
    end
    if ~strcmpi(allUnitNames2{unit},allUnitNames{unit})
        error('nont matching names')
    end
end

%%
%---Get Number of Parameters---%
numBinParams = [];
if any(contains(mCCVParms.variablesNames,'P'))
    numBinParams = [numBinParams binSizeParams.numSpatialBins^2];
end
if any(contains(mCCVParms.variablesNames,'D'))
    numBinParams = [numBinParams binSizeParams.numDirBins];
end
if any(contains(mCCVParms.variablesNames,'T'))
    numBinParams = [numBinParams binSizeParams.numTimesBins];
end
if any(contains(mCCVParms.variablesNames,'N'))
    numBinParams = [numBinParams 2];
end
if any(contains(mCCVParms.variablesNames,'E'))
    numBinParams = [numBinParams binSizeParams.numTimesBins]; %could have less coverage...
end

%%
%---Determine which Models Correspond to Which Parameters---%
[A,modelTypes,modelNames,numModels,numModelParameters] = ...
    buildModelDataMonkeyHpc(cell(1,length(mCCVParms.variablesNames)),mCCVParms.variablesNames);
modelTypes2 = cell2mat(modelTypes);
spatialModelsNumbers = find(modelTypes2(:,1) == 1);
directionModelsNumbers = find(modelTypes2(:,2) == 1);

%%
%---Get Time course of Beta Weights---%
numUnits = length(allMonkeyFixationLockedFiringSmoothed);
allPredictedRateMaps = cell(1,numUnits);
allPredictedResponses = cell(1,numUnits);
for unit  = 1:numUnits
    disp(['Processing Unit#' num2str(unit)])
    if isempty(allMeanParams{unit}) %no selectivity
        continue
    end

    %get linear-nonlinear poisson model filters
    paramFullModel = allMeanParams{unit};
    posParam = paramFullModel(1:numBinParams);
    sacDirParam = paramFullModel(numBinParams(1)+1:numBinParams(1)+numBinParams(2));
    timeParam = paramFullModel(sum(numBinParams(1:2))+1:end);

    %get this channels table
    thisChannelTable = allMonkeyFixationDataSmoothed{unit};


    %---Reconstruct the Firing Rate Map---%
    %get firing rate position bins and direction bins
    [~,~,~,allSpatialBinInd,spacingX,spacingY] = binSpaceGlm(thisChannelTable.posX,thisChannelTable.posY,thisChannelTable.meanResponse,binSizeParams);
    [dirVec,~,~,dirBinInd] = binDirectionGlm(thisChannelTable.saccadeDir,thisChannelTable.meanResponse,binSizeParams);
    meanTimeBeta = mean(timeParam);

    %get reconstructed responses
    allReconstructedResponses = NaN(size(thisChannelTable,1),1);
    for fix = 1:size(thisChannelTable,1)
        allReconstructedResponses(fix) = exp(posParam(allSpatialBinInd(fix)) + ...
            sacDirParam(dirBinInd(fix)) + meanTimeBeta);
    end
    %rescale since we don't know what the offset parameter is in the model
    allReconstructedResponses = allReconstructedResponses*mean(thisChannelTable.meanResponse)/mean(allReconstructedResponses);

    %get reconstructed spatial map aka predicted map
    reconstructedSpatialMap = NaN(1,numBinParams(1));
    for bin = 1:numBinParams(1)
        reconstructedSpatialMap(bin) = mean(allReconstructedResponses(allSpatialBinInd == bin));
    end
    reconstructedSpatialMap2 = reshape(reconstructedSpatialMap,sqrt(numBinParams(1)),sqrt(numBinParams(1)));


    %get predicted spatial
    predictedSpatialResponse = NaN(numBinParams(1),numBinParams(2),numBinParams(3));
    for binSpatial = 1:numBinParams(1)
        for binDir = 1:numBinParams(2)
            for binTime = 1:numBinParams(3)
                predictedSpatialResponse(binSpatial,binDir,binTime) =  exp(posParam(binSpatial) + ....
                    sacDirParam(binDir) + timeParam(binTime));
            end
        end
    end
    allPredictedResponses{unit} = predictedSpatialResponse;
    predictedSpatialResponse = mean(mean(predictedSpatialResponse,3),2);
    predictedSpatialResponse = predictedSpatialResponse*mean(thisChannelTable.meanResponse)/mean(predictedSpatialResponse);
    predictedSpatialResponse = reshape(predictedSpatialResponse,sqrt(numBinParams(1)),sqrt(numBinParams(1)));
    allPredictedRateMaps{unit} = predictedSpatialResponse;

    %---Plot Results for this Unit---%
    if plotFlag
        %%
        figure
        vals = allFiringRateMapsSmoothed{unit};
        subplot(2,4,1)
        imagesc(allFiringRateMapsSmoothed{unit},'AlphaData',~isnan(allFiringRateMapsSmoothed{unit}))
        axis off
        title('Spatial Tuning Curve')
        caxis([prctile(allFiringRateMapsSmoothed{unit}(:),2.5) prctile( allFiringRateMapsSmoothed{unit}(:),97.5)])
        colorbar

        subplot(2,4,2)
        imagesc(reshape(exp(posParam),sqrt(numBinParams(1)),sqrt(numBinParams(1))));
        axis off
        title('LNL Fit Response (exp(Beta))')
        colorbar

        subplot(2,4,3)
        imagesc(reconstructedSpatialMap2,'AlphaData',~isnan(reconstructedSpatialMap2))
        axis off
        title('Reconstructed Spatial Tuning Curve')
        colorbar

        subplot(2,4,4)
        imagesc(predictedSpatialResponse)
        axis off
        colorbar
        title('Predicted Spatial Map')

        subplot(2,4,5)
        plot(-40:40:400,exp(timeParam)/max(exp(timeParam)))
        hold on
        meanTime = mean(allMonkeyFixationLockedFiringSmoothed{unit});
        meanTime = meanTime/max(meanTime);
        plot(erpParameters.timeWindow,meanTime)
        hold off
        xlabel('Time from Fixation Onset')
        title('Scaled Temporal Responses')

        subplot(2,4,6)
        plot(dirVec*180/pi,exp(sacDirParam))
        title('Model Fit Saccade Direction')

        sgtitle(allUnitNames{unit})
    end
end


%%
%---Identify spatial model types--%
isSpatialModel = zeros(1,length(allBestModels));
for sMN = 1:length(spatialModelsNumbers)
    isSpatialModel(spatialModelsNumbers(sMN) == allBestModels) = 1;
end

%%
%---Combine All Models Together---%
allSpatialPredictedMaps = [];
allSpatialResponses = [];
allNonSpatialPredictedMaps = [];
allNonSpatialResponses = [];
for unit = 1:length(isSpatialModel)
    if isempty(allMeanParams{unit})
        continue
    end
    if isSpatialModel(unit) == 1
        allSpatialPredictedMaps = cat(3,allSpatialPredictedMaps, allPredictedRateMaps{unit});
        allSpatialResponses = cat(3,allSpatialResponses,squeeze(mean(allPredictedResponses{unit},2)));
    else
        allNonSpatialPredictedMaps = cat(3,allNonSpatialPredictedMaps, allPredictedRateMaps{unit});
        allNonSpatialResponses = cat(3,allNonSpatialResponses,squeeze(mean(allPredictedResponses{unit},2)));
    end
end


%% Population Decoding from Model Responses
%---Full Population Decoding---%
%for spatial units
spatialMlVal = NaN(numBinParams(1),numBinParams(3));
spatialBinErr = NaN(numBinParams(1),numBinParams(3));
spatialPostAvg = zeros(sqrt(numBinParams(1)),sqrt(numBinParams(1)),numBinParams(1));
for bin = 1:numBinParams(1)
    [ post ] = decode_calcPost(squeeze(allSpatialResponses(bin,:,:))', allSpatialPredictedMaps, tBin,true);
    [spatialBinErr(bin,:),spatialMlVal(bin,:)] = decodePost(post,bin);
    spatialPostAvg(:,:,bin) = spatialPostAvg(:,:,bin) + mean(post,3);
end


%---Spatial Limited Population---%
numNonSpatial = size(allNonSpatialResponses,3);
medianDecodingErrorLimited = NaN(1,numShuffs);
temporalDecodingError = NaN(numShuffs,numBinParams(3));
spatialDecodingErrorLimited = zeros(numBinParams(1),1);
for shuff = 1:numShuffs
    selectedChannels = randperm(size(allSpatialResponses,3));
    selectedChannels = selectedChannels(1:numNonSpatial);
    spatialMlValLimited = NaN(numBinParams(1),numBinParams(3));
    spatialBinErrLimited = NaN(numBinParams(1),numBinParams(3));
    for bin = 1:numBinParams(1)
        theseData = squeeze(allSpatialResponses(bin,:,selectedChannels))';
        [ post ] = decode_calcPost(theseData, allSpatialPredictedMaps(:,:,selectedChannels), tBin,true);
        [spatialBinErrLimited(bin,:),spatialMlValLimited(bin,:)] = decodePost(post,bin);
    end
    medianDecodingErrorLimited(shuff) = median(spatialBinErrLimited(:));
    temporalDecodingError(shuff,:) = mean(spatialBinErrLimited,1);
    spatialDecodingErrorLimited = spatialDecodingErrorLimited + mean(spatialBinErrLimited,2);
end
spatialDecodingErrorLimited = spatialDecodingErrorLimited/numShuffs;


%---For Full non-spatial units---%
nonSpatialMlVal = NaN(numBinParams(1),numBinParams(3));
nonSpatialBinErr = NaN(numBinParams(1),numBinParams(3));
nonSpatialPostAvg = zeros(sqrt(numBinParams(1)),sqrt(numBinParams(1)),numBinParams(1));
count = 0;
for bin = 1:numBinParams(1)
    [ post ] = decode_calcPost(squeeze(allNonSpatialResponses(bin,:,:))', allNonSpatialPredictedMaps,tBin,true);
    [nonSpatialBinErr(bin,:),nonSpatialMlVal(bin,:)] = decodePost(post,bin);
    nonSpatialPostAvg = nonSpatialPostAvg + mean(post,3);
end
medianNonSpatialDecodingError = median(nonSpatialBinErr(:));


%%
%---Run Stats---%
%Run 1-way ANOVA over Time
%data are close enough to normal because of meaning across spatial bins
[pAnova,tAnova,anovaStats] = anova1(temporalDecodingError);
multCompareResults = multcompare(anovaStats);

%run permutation tests on overal decoding accuracy
pDecodingAccuracy = ranksum(spatialBinErr(:),nonSpatialBinErr(:));
pShuffledAccuracy = 1-sum(medianDecodingErrorLimited< medianNonSpatialDecodingError)/numShuffs;

%%
%---Plot Decoding Results---%
%13.951 is the 97.5% decodeing error for the exact non-spatial units
figure
subplot(2,3,1)
imagesc(reshape(mean(spatialBinErr,2),sqrt(numBinParams(1)),sqrt(numBinParams(1))))
caxis([0 14])
axis off
title('Exact Spatial Neurons Decoding Error')

subplot(2,3,2)
imagesc(reshape(mean(nonSpatialBinErr,2),sqrt(numBinParams(1)),sqrt(numBinParams(1))))
caxis([0 14])
colorbar
axis off
title('Exact Non-Spatial Neurons Decoding Error')



subplot(2,3,3)
histogram(spatialBinErr(:),0:25)
hold on
histogram(nonSpatialBinErr(:),0:25)
hold off
box off
title(['Exact Decoding Error Histogram (p_{ranksum} = ' num2str(pDecodingAccuracy) ', p_{shuffled} = ' num2str(pShuffledAccuracy) ')'])
legend('Spatial','Non-Spatial')
ylabel('Bin Counts')

subplot(2,3,4)
imagesc(reshape(spatialDecodingErrorLimited,sqrt(numBinParams(1)),sqrt(numBinParams(1))))
axis off
colorbar
caxis([0 14])
title('Spatial Neurons Shuffled Average')

subplot(2,3,5)
hold on
plot(mean(temporalDecodingError),'k')
errorbar(mean(temporalDecodingError),std(temporalDecodingError)./sqrt(numShuffs),'k')
plot(mean(spatialBinErr,1))
hold off
legend({'Shuffled','','Exact'})
xlabel('Time Bin')
ylabel('Decoding Error')
title('Spatial: Shuffled Decoding Error Over Eye Time')
box off

subplot(2,3,6)
plot(mean(spatialBinErr,1))
hold on
plot(mean(nonSpatialBinErr,1))
hold off
xlabel('Time Bin')
ylabel('Decoding Accuracy')
legend('Spatial','Non-Spatial')
title('Exact Decoding Error Over Eye Time')
box off

sgtitle('Model Population-level Decoding')
save_and_close_fig(mainDataDir,'populationDecoding')
%%
%---Plot Decoding Over Time---%
figure
for i = 1:numBinParams(3)
    subplot(3,4,i)
    imagesc(reshape(spatialBinErr(:,i),sqrt(numBinParams(1)),sqrt(numBinParams(1))))
    title(mean(spatialBinErr(:,i)))
    axis off
    caxis([0 30])
end
sgtitle('Exact Spatial: Decoding Error by Time Window')
save_and_close_fig(mainDataDir,'populationspatialExactTimeCourse')

figure
for i = 1:numBinParams(3)
    subplot(3,4,i)
    imagesc(reshape(nonSpatialBinErr(:,i),sqrt(numBinParams(1)),sqrt(numBinParams(1))))
    title(mean(nonSpatialBinErr(:,i)))
    axis off
    caxis([0 30])
end
sgtitle('Exact Non-Spatial: Decoding Error by Time Window')
save_and_close_fig(mainDataDir,'populationNonSpatialExactTimeCourse')

%%
%---Save the Data---%
save([mainDataDir 'populationDecoding.mat'],'erpParameters','numShuffs','numBinParams','tBin',...
    'allFiringRateMapsSmoothed','allMonkeyFixationDataSmoothed','allMonkeyFixationLockedFiringSmoothed','allUnitNames',...
    'allSpatialPredictedMaps','allSpatialResponses','allNonSpatialPredictedMaps','allNonSpatialResponses',...
    'spatialMlVal','spatialBinErr','spatialPostAvg','nonSpatialMlVal','nonSpatialBinErr','nonSpatialPostAvg','medianNonSpatialDecodingError',...
    'numNonSpatial','medianDecodingErrorLimited','temporalDecodingError','spatialDecodingErrorLimited',...
    'pAnova','tAnova','anovaStats','multCompareResults',...
    'pDecodingAccuracy','pShuffledAccuracy');

%%
toc