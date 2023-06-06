function [A,modelTypes,modelNames,numModels,numModelParameters] = ...
    buildModelDataMonkeyHpc(dataCell,variableNames)
%function flexibly builds all model matrices of all combinations.
%Written by Seth Konig 7/25/22

%---Process Inputs---%
%get number of data types
numDataTypes = length(dataCell);

%determine number of models, sorted from most complicated to least
modelCombos = [];
for nDT = numDataTypes:-1:1
    theseCombos = nchoosek(1:numDataTypes,nDT);
    theseCombos2 = zeros(size(theseCombos,1),numDataTypes);
    for tC = 1:size(theseCombos,1)
        for tCC = 1:size(theseCombos,2)
            theseCombos2(tC,theseCombos(tC,tCC)) = 1;
        end
    end
    modelCombos = [modelCombos; theseCombos2];
end

%get some info about these models
numModelParameters = sum(modelCombos,2);%get number of data types in each model
numModels = size(modelCombos,1);


%---Build Data Matrices---%
A = cell(numModels,1); %all the data together
modelTypes = cell(numModels,1);
modelNames = cell(numModels,1);
for nM = 1:numModels
    modelTypes{nM} = modelCombos(nM,:);
    dataVariableNumbers = find(modelCombos(nM,:) == 1);
    A{nM} = [];
    modelNames{nM} = [];
    for dVN = 1:length(dataVariableNumbers)
         A{nM} = [A{nM} dataCell{dataVariableNumbers(dVN)}];
         modelNames{nM} = [modelNames{nM} variableNames{dataVariableNumbers(dVN)}];
    end
end

end