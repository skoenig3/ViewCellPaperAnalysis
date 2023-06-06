function [thisRandomNames,thisRandomEffectType,thisRandomEffects] = getRandomEffectStats(glme)
%own subfunction to keep it clean
%if more than random effect then includes corr

%get random effects
[~,~,randomStats] = covarianceParameters(glme);
if size(randomStats,1) == 1 %then error term only
    thisRandomNames = [];
    thisRandomEffectType = [];
    thisRandomEffects = [];
    return
end
numRandomEffects = size(randomStats,1)-1;

%store in more convient format
thisRandomEffectType = cell(1,numRandomEffects);
thisRandomNames = cell(1,numRandomEffects);
thisRandomEffects = NaN(1,numRandomEffects);
for nRE = 1:numRandomEffects
    thisRandomEffectType{nRE} = randomStats{nRE}.Type{1};
    thisRandomNames{nRE} = [randomStats{nRE}.Group ': ' randomStats{nRE}.Name1{1} ' vs ' randomStats{nRE}.Name2{1}];
    thisRandomEffects(nRE) =  randomStats{nRE}.Estimate(1);
end

end