% getCatStepVals.m
%
% Helper function for saveLegStepParamByCond_fly()
% Returns step values, as well as leg calls for selected amp and
%  dur, given keys. Operates on single vectors
%
% INPUTS:
%   stepValAll - vector of all step param vals, length n
%   whichAmp - amplitude value to select for
%   whichDur - duration value to select for
%   stepWhichLeg - vector of leg indices for all steps, length n
%   legStepCat - category each step belongs to, as indices into keys,
%     length n; m different indices
%   ampsKey - amplitude corresponding to each index in legStepCat, length m
%   dursKey - duration corresponding to each index in legStepCat, length m
%
% OUTPUTS:
%   stepVal - vector of step X positions, for this amp and dur, length l
%   theseWhichLeg - vector of leg indices, for this amp and dur, length l
%
% CREATED: 2/14/23 - HHY
%
% UPDATED:
%   2/14/23 - HHY
%
function [stepVal, theseWhichLeg] = getCatStepVals(stepValAll, ...
    whichAmp, whichDur, stepWhichLeg, legStepCat, ampsKey, dursKey)

    % turn whichAmp and whichDur into index value for legStepCat
    thisInd = intersect(find(ampsKey == whichAmp), ...
        find(dursKey == whichDur));

    % get indices into steps that match this categories
    valInd = find(legStepCat == thisInd);

    % use indices to select values of extreme point, which leg that are
    %  valid
    stepVal = stepValAll(valInd);
    theseWhichLeg = stepWhichLeg(valInd);   
end