function [timeVec,timeGrid,meanTimeResponse,timeCoverage] = binTimeGlm(timeValues,response,binSizeParams)
%function bins time values for glme
%could be used for other vector representations too!
%written by Seth Konig 7/25/2022

timeStep = round(max(timeValues)/binSizeParams.numTimesBins); %splits based on what we have
timeVec = timeStep:timeStep:max(timeValues);
newNumTimeBins = length(timeVec);
timeCoverage = zeros(1,newNumTimeBins);
timeGrid = zeros(size(timeValues,1),newNumTimeBins);
meanTimeResponse = zeros(1,newNumTimeBins);
for idx = 1:size(timeValues,1)
    [~, timeBin] = min(abs(timeValues(idx)-timeVec));
    timeGrid(idx, timeBin) = 1;
    timeCoverage(timeBin) =  timeCoverage(timeBin)+1;
    meanTimeResponse(timeBin) =  meanTimeResponse(timeBin) + response(idx);
end
meanTimeResponse = meanTimeResponse./timeCoverage;
end