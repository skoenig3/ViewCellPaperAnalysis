function [testFit,trainFit,paramMat] = fitGlmeMCCV(A,spiketrain,modelType,mCCVParms,numParams)
%Function is based off of fit_model_kfold from  Hardcastle, 2017 but
%uses Monte-Carlo Cross-Validation method instead of k-fold
%cross-validation because we have less data and there is biased data
%structure. Also we randomly seed every time.
%Written by Seth Konig 7/25/2022
%
% The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.


%---Initialize Everything---%
%get numbers and sizes
numDataPoints = size(A,1);
numCol = size(A,2);
numTrainInd = round(numDataPoints*mCCVParms.dataSplit);%e.g. 70%



%---Perform Monte-Carlo Cross Validation---%
testFit = nan(mCCVParms.shuffs,2);
trainFit = nan(mCCVParms.shuffs,2); 
paramMat = nan(mCCVParms.shuffs,numCol);
for shuffs = 1:mCCVParms.shuffs
    %fprintf('\t\t- Cross validation fold %d of %d\n', shuffs, mCCVParms.shuffs);
    
    %randomly assign data into testing and training proportions
    allInd = randperm(numDataPoints)';
    train_ind =  allInd(1:numTrainInd);
    test_ind  = allInd(numTrainInd+1:end);

    %get test data
    test_spikes = spiketrain(test_ind); %test spiking
    test_A = A(test_ind,:);
    
    % get training data
    train_spikes = spiketrain(train_ind);
    train_A = A(train_ind,:);
    
    %get initial parameters
    data{1} = train_A;
    data{2} = train_spikes;
    initParam = 1e-3*randn(numCol, 1); %randomly seed every time
    
    %fit model
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','off');
    param = fminunc(@(param) lnPoissonGlm(param,data,modelType,mCCVParms,numParams),initParam,options);
    
    % save the parameters
    paramMat(shuffs,:) = param;
    
    %Get Model Fitness
    [testFit(shuffs,1),testFit(shuffs,2)] = compute_llh(test_A,test_spikes,param);
    [trainFit(shuffs,1), trainFit(shuffs,2)] = compute_llh(train_A,train_spikes,param);

end

end