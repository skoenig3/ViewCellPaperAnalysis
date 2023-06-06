function [selectedModelForward,llhChangeValues,llhValues] = determineBestModelGlm(testFit,modelTypes,numDataPoints)
% Original code implemented forward feature selection in order to determine
% the simplest model that best describes neural spiking. First, the
% highest-performing single-variable model is identified. Then, the
% highest-perfmoring double-variable model that includes the
% single-variable model is identified. This continues until the full model
% is identified. Next, statistical tests are applied to see if including
% extra variables significantly improves model performance. The first time
% that including variable does NOT signficantly improve performance, the
% procedure is stopped and the model at that point is recorded as the
% selected model.


%---Identify Number of Levels---%
modelTypes2 = cell2mat(modelTypes);
numVariables = sum(modelTypes2,2);
maxLevel = max(numVariables);
if maxLevel > 4
    error('need code for this!')
end


%---Get Logliklihood Values---%
numShuffs = size(testFit{1},1);
numModels = size(testFit,1);
llhChangeValues = NaN(numShuffs,numModels); %what was LLH_values
llhValues = NaN(numShuffs,numModels);  %absolute LLH values
for nM = 1:numModels
    llhChangeValues(:,nM) = testFit{nM}(:,1);
    llhValues(:,nM) = testFit{nM}(:,2);
end

%get means
meanllhChangeValues = mean(llhChangeValues);
meanllhValues = mean(llhValues);



%---Forward Search---%
% find the best single model
singleModels = find(numVariables == 1);
[~,top1] = max(meanllhChangeValues(singleModels));
top1 = top1 + singleModels(1)-1;
bestVariable = find(modelTypes2(top1,:) == 1);

% find the best double model that includes the single model
doubleModels = find(numVariables == 2 & modelTypes2(:,bestVariable) == 1);
[~,top2] = max(meanllhChangeValues(doubleModels));
top2 = top2 + doubleModels(1)-1;
bestVariable2 = find(modelTypes2(top2,:) == 1);

% find the best triple model that includes the double model
if maxLevel >= 3
    trippleModels = find(numVariables == 3 & all(modelTypes2(:,bestVariable2) == 1,2));
    [~,top3] = max(meanllhChangeValues(trippleModels));
    top3 = top3 + trippleModels(1)-1;
end

%by default top 4 is the best top4
if maxLevel >= 4
    top4 = 1;
end


%sign-rank test to determine if models are better or not
if maxLevel == 4
    [pLlh12,~] = signrank(llhChangeValues(:,top2),llhChangeValues(:,top1),'tail','right');
    [pLlh23,~] = signrank(llhChangeValues(:,top3),llhChangeValues(:,top2),'tail','right');
    [pLlh34,~] = signrank(llhChangeValues(:,top4),llhChangeValues(:,top3),'tail','right');
elseif maxLevel == 3
    [pLlh12,~] = signrank(llhChangeValues(:,top2),llhChangeValues(:,top1),'tail','right');
    [pLlh23,~] = signrank(llhChangeValues(:,top3),llhChangeValues(:,top2),'tail','right');
    pLlh34 = 1;
elseif maxLevel == 2
    [pLlh12,~] = signrank(llhChangeValues(:,top2),llhChangeValues(:,top1),'tail','right');
    pLlh23 = 1;
    pLlh34 = 1;
    
end

%select the best model
if pLlh12 < 0.05 % double model is sig. better
    if pLlh23 < 0.05  % triple model is sig. better
        if pLlh34 < 0.05 % full model is sig. better
            selectedModelForward = 1; % full model
        else
            selectedModelForward = top3; %triple model
        end
    else
        selectedModelForward = top2; %double model
    end
else
    selectedModelForward = top1; %single model
end

% re-set if selected model is not above baseline
pValBaseline = signrank(llhChangeValues(:,selectedModelForward),[],'tail','right');
if pValBaseline > 0.05
    selectedModelForward = -1;
end


end