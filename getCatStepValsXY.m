% getCatStepValsXY.m
%
% Helper function for plotAEPPEPscatter(), plotAEPPEPContours(), and
%  plotAEPPEPmeans()
% Returns step X and Y values, as well as leg calls for selected amp and
%  dur, given keys. Operates on single vectors
%
% INPUTS:
%   stepEPX - vector of all step X extreme positions, length n
%   stepEPY - vector of all step Y extreme positions, length n
%   whichAmp - amplitude value to select for
%   whichDur - duration value to select for
%   stepWhichLeg - vector of leg indices for all steps, length n
%   legStepCat - category each step belongs to, as indices into keys,
%     length n; m different indices
%   ampsKey - amplitude corresponding to each index in legStepCat, length m
%   dursKey - duration corresponding to each index in legStepCat, length m
%
% OUTPUTS:
%   stepValX - vector of step X positions, for this amp and dur, length l
%   stepValY - vector of step Y positions, for this amp and dur, length l
%   theseWhichLeg - vector of leg indices, for this amp and dur, length l
%
% CREATED: 10/4/22 - HHY
%
% UPDATED:
%   10/4/22 - HHY
%
function [stepValX, stepValY, theseWhichLeg] = getCatStepValsXY(stepEPX, ...
    stepEPY, whichAmp, whichDur, stepWhichLeg, legStepCat, ampsKey, ...
    dursKey)

    % turn whichAmp and whichDur into index value for legStepCat
    thisInd = intersect(find(ampsKey == whichAmp), ...
        find(dursKey == whichDur));

    % get indices into steps that match this categories
    valInd = find(legStepCat == thisInd);

    % use indices to select values of extreme point, which leg that are
    %  valid
    stepValX = stepEPX(valInd);
    stepValY = stepEPY(valInd);
    theseWhichLeg = stepWhichLeg(valInd);   
end