% rmRefrPrdThreshCross.m
%
% Helper function called by detectSpikes() that removes threshold crossings
%  that occur within refractory period.
%
% INPUTS:
%   threshCrossInd  - indicies of threshold crossings
%   refrPrd - refractory period length, in samples
%
% OUTPUTS:
%   corrThreshCrossInd - indicies of threshold crossings, with invalid ones
%       removed
%
% CREATED: 9/13/20 - HHY
%
% UPDATED:
%   9/13/20 - HHY
%
function corrThreshCrossInd = rmRefrPrdThreshCross(threshCrossInd, refrPrd)

    % difference between each threshold crossing, in samples
    diffThreshInd = diff(threshCrossInd);
    
    % indicies of extra threshold crossings 
    extraCrossings = find(diffThreshInd < refrPrd) + 1;
    
    % remove extra threshold crossings (assumes one that follows too
    %  closely is wrong)
    corrThreshCrossInd = threshCrossInd;
    corrThreshCrossInd(extraCrossings) = [];
end