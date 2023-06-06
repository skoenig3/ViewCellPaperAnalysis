function [binErr,mlVal] = decodePost(post,allSpatialBinInd)
%based on decode_processPost but easier and faster
%Use ML  to decode posterior prob

%---Decode Bin Position---%
%select max bin ind
nTwin = size(post,3);
decodedRow = NaN(1,nTwin);
decodedCol = NaN(1,nTwin);
mlVal = NaN(1,nTwin);
for n = 1:nTwin
    curPost = post(:,:,n);
    [maxRow,maxCol ] = find(max(curPost(:)) == curPost);
    tmpInd = 1:length(maxRow);
    tmpInd = tmpInd(randperm(length(tmpInd),1));%If multiple take 1@rand
    maxRow = maxRow(tmpInd);
    maxCol = maxCol(tmpInd);
    mlVal(n) = curPost(maxRow(1),maxCol(1)); %Store value in post at ML decode location
    decodedRow(n) = maxRow(1);
    decodedCol(n) = maxCol(1);
end

%---Get Bin Error---%
%get bin x/y
[trueY,trueX] = ind2sub([25 25],allSpatialBinInd);
binErr = sqrt((trueY-decodedRow).^2 + (trueX-decodedCol).^2);

end