% getCatStepValsCond.m
%
% Helper function for saveLegStepParamByCond_fly()
% Returns step val, as well as leg calls for selected amp and
%  dur, given keys. Conditions on step value. Operates on single vectors. 
% Modification of getCatStepValsXYcond.m
%
% INPUTS:
%   stepValAll - vector of all step param vals, length n
%   whichAmp - amplitude value to select for
%   whichDur - duration value to select for
%   condVec - vector of all step values to condition on, length n
%   condState - string describing statement to condition on
%   stepWhichLeg - vector of leg indices for all steps, length n
%   legStepCat - category each step belongs to, as indices into keys,
%     length n; m different indices
%   ampsKey - amplitude corresponding to each index in legStepCat, length m
%   dursKey - duration corresponding to each index in legStepCat, length m
%
% OUTPUTS:
%   stepVal - vector of step param vals, for this amp and dur, length l
%   theseWhichLeg - vector of leg indices, for this amp and dur, length l
%
% CREATED: 2/14/23 - HHY
%
% UPDATED:
%   2/14/23 - HHY
%
function [stepVal, theseWhichLeg] = getCatStepValsCond(stepValAll, ...
    whichAmp, whichDur, condVec, condState, stepWhichLeg, ...
    legStepCat, ampsKey, dursKey)

    % turn whichAmp and whichDur into index value for legStepCat
    thisInd = intersect(find(ampsKey == whichAmp), ...
        find(dursKey == whichDur));

    % get indices into steps that match this categories
    valInd = find(legStepCat == thisInd);

    % use indices to select values of extreme point, which leg that are
    %  valid
    tempStepVal = stepValAll(valInd);
    tempTheseWhichLeg = stepWhichLeg(valInd);  
    tempCondVec = condVec(valInd);

    % of steps that meet amp and dur category, find ones that also meet
    %  condition
    finalValInd = find(eval(['tempCondVec' condState]));

    stepVal = tempStepVal(finalValInd);
    theseWhichLeg = tempTheseWhichLeg(finalValInd);
end