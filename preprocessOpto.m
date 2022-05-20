% preprocessOpto.m
%
% Function for processing raw voltage signals indicating when the
%  optogenetic stimulation was active (feedback signal from shutter in 
%  front of mercury lamp; high = open). Extracts start and end times of
%  opto stim. Also, logical for stim on or not. Saves parameters of opto
%  stim that were originally saved in inputParams. 
% Returns these parameters in opto struct
%
% INPUTS:
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqOutput - struct of output signals sent by experimental DAQ,
%       processed by preprocessUserDaq() 
%   daqTime - vector of times corresponding to each sample point of daqData
%   inputParams - struct of input parameters, specific to particular
%       experiment type, saved in trial.mat file when experiment run
%
% OUTPUTS:
%   opto - output struct with the following fields:
%       stimStartTimes - start times of the opto stim, in seconds
%       stimEndTimes - end times of the opto stim, in seconds
%       stimDurs - actual duration of each opto stim bout, in seconds
%       stimCmdDurs - commanded duration, in seconds (for each stim, nearest
%           commanded duration)
%       stimOnLogical - logical for whether opto stim is on (1) or not (0);
%           num elements from sample rate of DAQ
%       stimParams - struct of parameters of opto stim, from inputParams
%
% CREATED: 3/28/22 - HHY
%
% UPDATED:
%   3/28/22 - HHY
%
function opto = preprocessOpto(daqData, daqTime, inputParams)

    % stim start indices, from shutter sync out signal; rising edges
    stimStarts = find(diff(daqData.HgLampShutterSyncOut) > 0.1);

    % stim end indices, from shutter sync out signal; falling edges
    stimEnds = find(diff(daqData.HgLampShutterSyncOut) < -0.1);

    % check that stimStarts and stimEnds are same length (should be except
    % in edge case where stim was still on when trial ended); if stimEnds
    % shorter than stimStarts, add in end (last index)
    if (length(stimEnds) < length(stimStarts))
        stimEnds(length(stimEnds) + 1) = ...
            length(daqData.HgLampShutterSyncOut) - 1;
    elseif (length(stimStarts) < length(stimEnds))
        disp('Warning: fewer opto stim starts than ends');
    end

    % stim start and end times
    stimStartTimes = daqTime(stimStarts + 1);
    stimEndTimes = daqTime(stimEnds + 1);

    % get logical for opto stim on
    stimOnLogical = daqData.HgLampShutterSyncOut > 0.5;

    % stimulation durations
    stimDurs = stimEndTimes - stimStartTimes;

    % convert actual stimulation durations to commanded durations
    cmdDurs = inputParams.optoStimParams.allStimDurs;
    stimCmdDurs = zeros(size(stimDurs)); % preallocate

    for i = 1:length(stimDurs)
        thisStimDur = stimDurs(i);

        % distance between this stimulation duration and all commanded
        %  durations
        durDist = abs(cmdDurs - thisStimDur);
        % get the nearest commanded duration
        [~, nearInd] = min(durDist);
        thisCmdDur = cmdDurs(nearInd);

        % update vector 
        stimCmdDurs(i) = thisCmdDur;
    end

    % add to output struct
    opto.stimStartTimes = stimStartTimes;
    opto.stimEndTimes = stimEndTimes;
    opto.stimDurs = stimDurs;
    opto.stimCmdDurs = stimCmdDurs;
    opto.stimOnLogical = stimOnLogical;

    % stimulation parameters
    opto.stimParams = inputParams.optoStimParams;
    opto.stimParams.optoStimProtocol = inputParams.optoStimProtocol;


end

