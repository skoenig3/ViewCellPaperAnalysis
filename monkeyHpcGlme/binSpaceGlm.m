function [posVec,posGrid,meanPosResponse,allBinInd,spacingX,spacingY] = binSpaceGlm(posX,posY,response,binSizeParams)
%function bins space in to square bins
%written by Seth Konig 7/25/2022

spacingX = binSizeParams.dimX/binSizeParams.numSpatialBins;
spacingY = binSizeParams.dimY/binSizeParams.numSpatialBins;
xbin = spacingX:spacingX:binSizeParams.dimX;
ybin = spacingY:spacingY:binSizeParams.dimY;
posVec ={xbin,ybin};
posGrid = zeros(size(posX,1),binSizeParams.numSpatialBins.^2);
posCoverage = zeros(1,binSizeParams.numSpatialBins.^2);
meanPosResponse = zeros(1,binSizeParams.numSpatialBins.^2);
allBinInd = NaN(size(posX,1),1);
for idx = 1:size(posX,1)
    [~, xcoor] = min(abs(posX(idx)-xbin));
    [~, ycoor] = min(abs(posY(idx)-ybin));
    binIdx = sub2ind([binSizeParams.numSpatialBins, binSizeParams.numSpatialBins], binSizeParams.numSpatialBins - ycoor + 1, xcoor);
    posGrid(idx, binIdx) = 1;
    posCoverage(binIdx) = posCoverage(binIdx) + 1;
    allBinInd(idx) = binIdx;
    if ~isempty(response)
        meanPosResponse(binIdx) = meanPosResponse(binIdx) + response(idx);
    end
end
meanPosResponse = meanPosResponse./posCoverage;
end
