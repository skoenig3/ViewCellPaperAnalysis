function [dirVec,dirGrid,meanDirResponse,allBinIds] = binDirectionGlm(direction,response,binSizeParams)
%function bins directional values for glme
%written by Seth Konig 7/25/2022

dirVec = 2*pi/binSizeParams.numDirBins/2:2*pi/binSizeParams.numDirBins:2*pi-2*pi/binSizeParams.numDirBins/2;
dirGrid = zeros(size(direction,1),binSizeParams.numDirBins);
dirCoverage = zeros(1,binSizeParams.numDirBins);
meanDirResponse = zeros(1,binSizeParams.numDirBins);
allBinIds = NaN(size(direction,1),1);
for idx = 1:size(direction,1)
    [~, dirBin] = min(abs(direction(idx)-dirVec));
    dirGrid(idx, dirBin) = 1;
    dirCoverage(dirBin) =  dirCoverage(dirBin)+1;
    allBinIds(idx) = dirBin;
    if ~isempty(response)
        meanDirResponse(dirBin) = meanDirResponse(dirBin) + response(idx);
    end
end
meanDirResponse = meanDirResponse./dirCoverage;
end