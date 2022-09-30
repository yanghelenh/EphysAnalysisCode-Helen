% shiftIinjStartTime.m
%
% Function that takes in iInj struct (output of preprocessIInj() and an
%  amount to shift and returns a modified iInj struct with current
%  injection times shifted by input amount (and other struct fields
%  modified to match)
% Helper function for computeLegStepsFictracEphysIinj_1Fly()
%
% INPUTS:
%   iInj - iInj struct, output from preprocessIInj()
%   shiftTime - time to shift current injection start times, in sec;
%       postive shifts move start time later, negative shifts earlier
%
% OUTPUTS:
%   iInjMod - modified iInj struct, with following fields (subset of those
%       in original iInj, only those required by
%       computeLegStepsFictracEphysIinj_1Fly()):
%       startTimes - current injection start times, in sec, shifted by
%           shiftTime
%       endTimes - current injection end times, in sec, shifted by
%           shiftTime
%       amps - amplitude of each current injection step, num ele matches
%           that of startTimes (which might not be same as original)
%       durs - duration of each current injection step
%
% CREATED: 9/21/22 - HHY
%
% UPDATED:
%   9/22/22 - HHY
%
function iInjMod = shiftIinjStartTime(iInj, shiftTime)
    % get shifted startTimes, endTimes
    newStartTimes = iInj.startTimes + shiftTime;
    newEndTimes = iInj.endTimes + shiftTime;

    % create iInjMod struct with amps and durs fields
    iInjMod.amps = iInj.amps;
    iInjMod.durs = iInj.durs;

    % remove any negative startTimes (all other invalid start times get
    %  taken care of in extractTrialsIinj())
    invalidInd = find(newStartTimes < 0);

    % remove trials with invalid (negative) startTimes, for all vars
    if ~isempty(invalidInd)
        newStartTimes(invalidInd) = [];
        newEndTimes(invalidInd) = [];
        iInjMod.amps(invalidInd) = [];
        iInjMod.durs(invalidInd) = [];
    end

    % add start and end times to iInjMod
    iInjMod.startTimes = newStartTimes;
    iInjMod.endTimes = newEndTimes;


end