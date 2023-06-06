function [bestModels,allMeanParams] = runLNGlmeMonkeyHPC(figureDir,dataName,allFiringRateMaps,...
    allFixationData,mCCVParms,binSizeParams)
%written by Seth Konig 7/24/2022

%---Run Analysis On One Unit at a Time---%
allMeanParams = cell(1,length(allFixationData));
bestModels = NaN(1,length(allFixationData));
for unit = 1:length(allFixationData)
    %check if there's enough data t run on this unit
    if isempty(allFixationData{unit})
        continue
    elseif sum(allFixationData{unit}.meanResponse > 0) < 0.01*size(allFixationData{unit},1)
        disp('Unit has too few responses/activity is too low so skipping...')
        bestModels(1,unit) = -99;
        continue
    end
    
    
    %---Grab Data---%
    %grab this data
    thisFixationData = allFixationData{unit};
    
    %remove fixations before 500 ms as we normally do. also makes fitting hard
    tooEarly = thisFixationData.fixationTime < 500;
    thisFixationData(tooEarly,:) = [];
    
    %get firing rate data
    firingRate = thisFixationData.meanResponse;
    
    
    %---Reformat and Prepare The Data---%
    %store all parameters in 1 cell
    allData = {};
    numBinParams = [];
    
    %prepare spatial data and turn into matrix
    if any(contains(mCCVParms.variablesNames,'P'))
        [posVec,posGrid,meanPosResponse] = binSpaceGlm(thisFixationData.posX,thisFixationData.posY,firingRate,binSizeParams);
        allData = [allData {posGrid}];
        numBinParams = [numBinParams binSizeParams.numSpatialBins^2];
    end
    
    %prepare saccade direction and turn into matrix too
    if any(contains(mCCVParms.variablesNames,'D'))
        thisFixationData.saccadeDir = thisFixationData.saccadeDir*pi/180+pi;%convert to radians 1st
        [dirVec,dirGrid,meanDirResponse] = binDirectionGlm(thisFixationData.saccadeDir,firingRate,binSizeParams);
          allData = [allData {dirGrid}];
           numBinParams = [numBinParams binSizeParams.numDirBins];
    end
    
    %prepare time vectors and turn into matrix too
    if any(contains(mCCVParms.variablesNames,'T'))
        [timeVec,timeGrid,meanTimeResponse,timeCoverage] = binTimeGlm(thisFixationData.fixationTime,firingRate,binSizeParams);
             allData = [allData {timeGrid}];
             numBinParams = [numBinParams binSizeParams.numTimesBins];
    end
    
    %get novelty vector
    if any(contains(mCCVParms.variablesNames,'N'))
        noveltyValues = [0 1];%0 for novel, 1 for repeat
        noveltyVal = thisFixationData.novelRepeatImage-1;%not really a grid but keeping name convention
        noveltyGrid = [noveltyVal == 0 noveltyVal == 1];
        meanNoveltyResponse = [mean(firingRate(thisFixationData.novelRepeatImage == 1)) ...
            mean(firingRate(thisFixationData.novelRepeatImage == 2))];
        allData = [allData {noveltyGrid}];
        numBinParams = [numBinParams 2];
    end
    
    %for eye time
    if any(contains(mCCVParms.variablesNames,'E'))
       [eyeVec,eyeGrid,meanEyeResponse,eyeCoverage] = binTimeGlm(thisFixationData.eyePhase,firingRate,binSizeParams);
        allData = [allData {eyeGrid}];
        numBinParams = [numBinParams length(eyeCoverage)]; %could have less coverage...
    end
         
    
    
    %---Do LN Model Fit---%
    %put into flexible structure that's nice and easy to use
    [A,modelTypes,modelNames,numModels,numModelParameters] = ...
        buildModelDataMonkeyHpc(allData,mCCVParms.variablesNames);
    
    %fit models
    spiketrain = round(firingRate);%round to make interger for Poisson model (not ideal/kosher)
    testFit = cell(numModels,1);
    trainFit = cell(numModels,1);
    param = cell(numModels,1);
    for nM = 1:numModels
        fprintf('\t- Fitting model %d of %d\nM', nM, numModels);
        [testFit{nM},trainFit{nM},param{nM}] = fitGlmeMCCV(A{nM},spiketrain,modelTypes{nM},mCCVParms,numBinParams);
    end
    
    
    %---Find the Best Model---%
    numDataPoints = size(thisFixationData,1);
    [selectedModelForward,llhChangeValues,llhValues] = determineBestModelGlm(testFit,modelTypes,numDataPoints);
    
    %get mean and std for plotting
    llhIncreaseMean = mean(llhChangeValues);
    llhIncreaseSem = std(llhChangeValues)/sqrt(mCCVParms.shuffs);
    
    
    
    %---Plot the Results---%
    figure
    subplot(3,3,1)
    imagesc(allFiringRateMaps{unit},'AlphaData',~isnan(allFiringRateMaps{unit}))
    c = caxis;
    caxis([c(1) prctile(allFiringRateMaps{unit}(:),97.5)])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('Observed: Rate Map')
    colorbar
    
    subplot(3,3,2)
    plot(dirVec,meanDirResponse)
    xlabel('Direction (rad)')
    ylabel('Firing Rate (Hz)')
    title('Observed: Binned Saccade Direction')
    
    subplot(3,3,3)
    plot(eyeVec,meanEyeResponse)
    xlabel('Time (bin/phase)')
    ylabel('Firing Rate (Hz)')
    title('Observed: Eye Phase/Time')
    
    % show parameters from the full model
    paramFullModel = mean(param{1},1);
    
    % pull out the parameter values
    posParam = paramFullModel(1:numBinParams);
    sacDirParam = paramFullModel(numBinParams(1)+1:numBinParams(1)+numBinParams(2));
    timeParam = paramFullModel(sum(numBinParams(1:2))+1:end);

    % since poisson model, beta weights are defaulted in log(rate ratio) so
    % convert to rate ratio
    posBeta = exp(posParam);
    sacDirBeta = exp(sacDirParam);
    eyePhaseBeta = exp(timeParam);
    
    % plot the model-derived response profiles
    subplot(3,3,4)
    fittedResponse = reshape(posBeta,binSizeParams.numSpatialBins,binSizeParams.numSpatialBins);
    imagesc(fittedResponse);
    c =caxis;
    caxis([c(1) prctile(fittedResponse(:),97.5)])
    axis off;
    colorbar
    title('LN Fit: Beta_{spatial}')
    
    subplot(3,3,5)
    plot(dirVec,sacDirBeta,'k')
    xlabel('Direction (rad)')
    title('LN Fit: Beta_{Saccade Direction}')
    
    subplot(3,3,6)
    plot(eyeVec,eyePhaseBeta,'k')
    xlabel('Time (bin/phase)')
    title('LN Fit: Beta_{Eye Phase}')
   
    subplot(3,3,7:9)
    errorbar(llhIncreaseMean,llhIncreaseSem,'ok','linewidth',3)
    hold on
    if selectedModelForward > 0
        plot(selectedModelForward,llhIncreaseMean(selectedModelForward),'.r','markersize',25)
    end
    plot(0.5:15.5,zeros(16,1),'--b','linewidth',2)
    hold off
    box off
    set(gca,'XLim',[0 length(llhIncreaseMean)+1]);
    set(gca,'XTick',1:length(llhIncreaseMean))
    set(gca,'XTickLabel',modelNames);
    if selectedModelForward > 0
        legend({'Model performance','Selected model','Baseline'},'Location','northeastoutside')
        title(['Best Model was ' modelNames{selectedModelForward}])
    else
        legend({'Model performance','Best AIC','Best BIC','Baseline'},'Location','northeastoutside')
        title('No good Model!')
    end
    
    
    %save and close figure as png and .fig
    save_and_close_fig(figureDir,[dataName '_' num2str(unit) '_LnModelFit'])

    
    %---Store Data Across Units---%
    %and prepare for export
    allMeanParams{unit} = paramFullModel;
    bestModels(unit) = selectedModelForward;
end

end