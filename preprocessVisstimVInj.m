% preprocessVisstimVInj.m
%
% Function for processing voltage signals indicating what visual stimulus
%  was presented on the G3 visual arena. Operates specifically on visual
%  stimuli controlled in closed-loop through voltage injection from the DAQ
% Reads in stimulus name and voltage injection protocol names to interpret
%  voltage signals to extract visual stimulus start and stop times and
%  velocities
% Currently, only interprets a small subset of voltage protocols
%  (spacedMultiRampVInj and multiStepVInj) used together to generate an
%  optomotor . Add more if additional visual
%  stimuli used
%
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
%   visstim - output struct with the following fields:
%       visstimParams - parameters of visual stimulus with following fields
%           visstimMode
%           patternName - name of visual stimulus pattern
%           patternIndex - index of visual stimulus pattern
%           vInjProtocolX - name of vInj protocol for X
%       
%
% CREATED: 1/16/24 - HHY
%   
% UPDATED:
%   1/16/24 - HHY
%   1/18/24 - HHY - currently written to be quite specific for data
%       gathered between 1/10/24 and 1/15/24
%   4/2/24 - HHY - updated so it analyzes optomotor experiments from
%       1/26/24 onwards
%
function visstim = preprocessVisstimVInj(daqData, daqTime, inputParams, ...
    daqOutput)

    % CONSTANTS
    V_RANGE = 10; % voltage goes from 0 to 10
    DEG_PER_PX = 360/96; % degrees per pixel of arena
%     % indices of regular grating patterns
%     REG_GRATING_PATTERN_IND = [9 20 21]; 
    % indices of new grating visual stimuli (from 1/26/24 onwards)
    REG_GRATING_PATTERN_IND_NEW = [22];

    % check visual stimulus mode
    % closed-loop, controlled with voltage injection from DAQ
    if (strcmp(inputParams.visstimMode, 'closedLoopXY'))
        % add visual stimulus mode to output struct
        visstim.visstimParams.visstimMode = inputParams.visstimMode;

        % add pattern name and index to output struct
        visstim.visstimParams.patternName = inputParams.patternName{1};
        visstim.visstimParams.patternIndex = inputParams.patternIndex;
        
        % load in pattern used
        patternPath = [patternsDir() filesep inputParams.patternName{1}];
        load(patternPath, 'pattern');

%% analysis of optomotor experiments from 1/26/24 onwards
    
        % process X and Y signals of spacedMultiRampVInj and multiStepVInj 
        %  on appropriate grating patterns (part of list in constants)
        if (strcmp(inputParams.vInjProtocolY, 'multiStepVInj') && ...
                any(inputParams.patternIndex == REG_GRATING_PATTERN_IND_NEW) && ...
                strcmp(inputParams.vInjProtocolX, 'spacedMultiRampVInj'))
            % add visual stimulus X protocol to output struct
            visstim.visstimParams.vInjProtocolX = ...
                inputParams.vInjProtocolX;
            % add visual stimulus Y protocol to output struct
            visstim.visstimParams.vInjProtocolY = ...
                inputParams.vInjProtocolY;

            % add info about commands to output struct
            visstim.visstimParams.yCmdAmps = ...
                inputParams.vInjParamsY.stepAmps;
            visstim.visstimParams.yCmdDurs = ...
                inputParams.vInjParamsY.stepDurs;

            % number of possible voltage values in Y (will be 4)
            % 1 - grating (static), 2 - grating (ramp), 3 - gray, 4 - off
            numSteps = pattern.y_num;

            % possible commanded output voltages
            cmdOutV = (1:numSteps) * (V_RANGE/numSteps);

            % thresholds given possible output voltages
            % set at midpoint between output voltages
            stepVThresh = (cmdOutV(1:(end-1)) + cmdOutV(2:end))/2;

            % interpret steps given pattern
            % regular square wave gratings
            %  1 - grating (static), 2 - grating (ramp), 3 - gray, 4 - off

            % get logical and start and end times for when grating was
            %  static
            staticLogical = (daqData.panelsDAC1Y' < stepVThresh(1));

            % start and end indices for grating static
            staticStarts = find(diff(staticLogical) > 0.1);
            staticEnds = find(diff(staticLogical) < -0.1);

            if (~isempty(staticStarts) && ~isempty(staticEnds))
                % check that end indices follow start indices at beginning and
                %  end of the trial, remove if not
                % beginning
                % if first start is after first end, remove first end
                if (staticStarts(1) > staticEnds(1))
                    staticEnds(1) = [];
                end
                % end
                % if last start is after last end, remove last start
                if (staticStarts(end) > staticEnds(end))
                    staticStarts(end) = [];
                end
            end

            % start and end times for grating displayed
            staticStartTimes = daqTime(staticStarts);
            staticEndTimes = daqTime(staticEnds);

            % get logical and start and end times for when grating moving
            rampLogical = (daqData.panelsDAC1Y' > stepVThresh(1)) & ...
                (daqData.panelsDAC1Y' < stepVThresh(2));

            % start and end indices for grating moving (ramp)
            rampStarts = find(diff(rampLogical) > 0.1);
            rampEnds = find(diff(rampLogical) < -0.1);

            if (~isempty(rampStarts) && ~isempty(rampEnds))
                % check that end indices follow start indices at beginning and
                %  end of the trial, remove if not
                % beginning
                % if first start is after first end, remove first end
                if (rampStarts(1) > rampEnds(1))
                    rampEnds(1) = [];
                end
                % end
                % if last start is after last end, remove last start
                if (rampStarts(end) > rampEnds(end))
                    rampStarts(end) = [];
                end
            end

            % start and end times for grating moving
            rampStartTimes = daqTime(rampStarts);
            rampEndTimes = daqTime(rampEnds);

            % get logical and start and end times for when gray
            %  displayed
            grayLogical = (daqData.panelsDAC1Y' > stepVThresh(2)) & ...
                (daqData.panelsDAC1Y' < stepVThresh(3));

            % start and end indices for grating displayed
            grayStarts = find(diff(grayLogical) > 0.1);
            grayEnds = find(diff(grayLogical) < -0.1);

            if (~isempty(grayStarts) && ~isempty(grayEnds))
                % check that end indices follow start indices at beginning and
                %  end of the trial, remove if not
                % beginning
                % if first start is after first end, remove first end
                if (grayStarts(1) > grayEnds(1))
                    grayEnds(1) = [];
                end
                % end
                % if last start is after last end, remove last start
                if (grayStarts(end) > grayEnds(end))
                    grayStarts(end) = [];
                end
            end

            % start and end times for gray
            grayStartTimes = daqTime(grayStarts);
            grayEndTimes = daqTime(grayEnds);


            % compare durations of each part of trial to what was commanded
            %  if duration is less than 0.9x or more than 1.1x of what was 
            %  expected, remove that start-end pair

            % static = 1
            staticCmdDur = visstim.visstimParams.yCmdDurs(1);
            staticActDurs = staticEndTimes - staticStartTimes;
            staticRmvLog = (staticActDurs < (0.9 * staticCmdDur)) | ...
                (staticActDurs > (1.1 * staticCmdDur));
            % remove any static portions that are too short
            staticStarts(staticRmvLog) = [];
            staticStartTimes(staticRmvLog) = [];
            staticEnds(staticRmvLog) = [];
            staticEndTimes(staticRmvLog) = [];
            staticActDurs(staticRmvLog) = [];

            % ramp = 2
            rampCmdDur = visstim.visstimParams.yCmdDurs(2);
            rampActDurs = rampEndTimes - rampStartTimes;
            rampRmvLog = (rampActDurs < (0.9 * rampCmdDur)) | ...
                (rampActDurs > (1.1 * rampCmdDur));
            % remove any ramp portions that are too short/long
            rampStarts(rampRmvLog) = [];
            rampStartTimes(rampRmvLog) = [];
            rampEnds(rampRmvLog) = [];
            rampEndTimes(rampRmvLog) = [];

            % gray = 3
            grayCmdDur = visstim.visstimParams.yCmdDurs(3);
            grayActDurs = grayEndTimes - grayStartTimes;
            grayRmvLog = (grayActDurs < (0.9 * grayCmdDur)) | ...
                (grayActDurs > (1.1 * grayCmdDur));
            % remove any gray portions that are too short/long
            grayStarts(grayRmvLog) = [];
            grayStartTimes(grayRmvLog) = [];
            grayEnds(grayRmvLog) = [];
            grayEndTimes(grayRmvLog) = [];
            grayActDurs(grayRmvLog) = [];



            % add to output
            visstim.staticStartTimes = staticStartTimes;
            visstim.staticEndTimes = staticEndTimes;
            visstim.staticOnLogical = staticLogical;
            visstim.staticDurs = staticActDurs;
            visstim.rampStartTimes = rampStartTimes;
            visstim.rampEndTimes = rampEndTimes;
            visstim.rampOnLogical = rampLogical;
            visstim.grayStartTimes = grayStartTimes;
            visstim.grayEndTimes = grayEndTimes;
            visstim.grayOnLogical = grayLogical;
            visstim.grayDurs = grayActDurs;


            % get ramp velocities, in degrees/sec, using Y channel for ramp
            %  start and end times/indices
            % volts per pixel
            vPerPx = V_RANGE / pattern.x_num;

            % ramp durations
            rampDurs = rampEndTimes - rampStartTimes;

            % preallocate
            rampVel = zeros(size(rampDurs));

            % loop through each ramp, get velocity
            for i = 1:length(rampDurs)
                % ramp start and end indices
                thisRampStart = rampStarts(i);
                thisRampEnd = rampEnds(i);

                % max voltage and index of max voltage
                [maxV, maxVInd] = max(...
                    daqData.panelsDAC0X(thisRampStart:thisRampEnd));
                % min voltage and index of min voltage
                [minV, minVInd] = min(...
                    daqData.panelsDAC0X(thisRampStart:thisRampEnd));

                % number of LEDs moved
                numLEDs = round((maxV - minV) / vPerPx);

                % speed of this ramp
                thisRampSpd = (numLEDs * DEG_PER_PX) / rampDurs(i);

                % get direction of movement
                % if index of max is less than the index of the min
                %  stimulus rotating CW -> (+)
                if (maxVInd < minVInd)
                    rampVel(i) = thisRampSpd;
                % if index of max is greater than the index of the min 
                %  stimulus rotating CCW -> (-)
                else
                    rampVel(i) = thisRampSpd * -1;
                end
            end

            % convert actual ramp durations to commanded durations
            cmdRampDurs = unique(inputParams.vInjParamsX.rampDur);

            rampsCmdDurs = zeros(size(rampDurs)); % preallocate
        
            for i = 1:length(rampDurs)
                thisRampDur = rampDurs(i);
        
                % distance between this stimulation duration and all 
                %  commanded durations
                durDist = abs(cmdRampDurs - thisRampDur);
                % get the nearest commanded duration
                [~, nearInd] = min(durDist);
                thisCmdDur = cmdRampDurs(nearInd);
        
                % update vector 
                rampsCmdDurs(i) = thisCmdDur;
            end

            % convert actual ramp velocities to commanded velocities
            cmdVel = (((inputParams.vInjParamsX.startV - ...
                inputParams.vInjParamsX.endV)/vPerPx) * DEG_PER_PX) ./ ...
                inputParams.vInjParamsX.rampDur;

            rampsCmdVels = zeros(size(rampDurs)); % preallocate
            for i = 1:length(rampDurs)
                thisRampVel = rampVel(i);
        
                % distance between this rampVel and all 
                %  commanded vel
                velDist = abs(cmdVel - thisRampVel);
                % get the nearest commanded velocity
                [~, nearInd] = min(velDist);
                thisCmdVel = cmdVel(nearInd);
        
                % update vector 
                rampsCmdVels(i) = thisCmdVel;
            end

            % add to output struct
            visstim.rampDurs = rampDurs;
            visstim.rampCmdDurs = rampsCmdDurs;
            visstim.rampVels = rampVel;
            visstim.rampCmdVels = rampsCmdVels;

        end

%% analysis of optomotor experiments from 1/10/24 to 1/15/24
%  uses spacedMultiRampVInj for X and multiStepVInj for Y
%         % process X channel signal - different depending on VInj protocol
%         % spacedMultiRampVInj
%         if (strcmp(inputParams.vInjProtocolX, 'spacedMultiRampVInj'))
%             % add visual stimulus X protocol to output struct
%             visstim.visstimParams.vInjProtocolX = ...
%                 inputParams.vInjProtocolX;
% 
%             % assuming values during space b/w ramps are more common than
%             %  any ramp values (should be the case unless space very short)
%             %  get voltage value during space
%             % assumes b/f and a/f space voltages are the same
%             spaceV = mode(daqData.panelsDAC0X);
%             % threshold for when ramp starts/ends
%             % assumes space voltages commanded at 0
%             rampVThresh = spaceV * 1.5;
% 
%             % logical for when voltage above threshold (during ramp)
%             rawRampLog = daqOutput.ampExtCmdIn > rampVThresh;
% 
%             % get ramp starts (rising edges in logical)
%             rampStarts = find(diff(rawRampLog) > 0.1);
% 
%             % get ramp ends (falling edges in logical)
%             rampEnds = find(diff(rawRampLog) < -0.1);
% 
%             % check that rampStarts and rampEnds are same length (should 
%             %  be except in edge case where ramp was still on when trial 
%             %  ended); if rampEnds shorter than rampStarts, add in end 
%             %  (last index)
%             if (length(rampEnds) < length(rampStarts))
%                 rampEnds(length(rampEnds) + 1) = ...
%                     length(daqOutput.ampExtCmdIn) - 1;
%             elseif (length(rampStarts) < length(rampEnds))
%                 disp('Warning: fewer ramp starts than ends');
%             end
% 
%             % convert indices to times
%             rampStartTimes = daqTime(rampStarts + 1);
%             rampEndTimes = daqTime(rampEnds + 1);
% 
%             rampDurs = rampEndTimes - rampStartTimes;
% 
%             % logical for when ramp is on
%             rampOnLogical = false(size(daqTime));
% 
%             % loop through all ramps
%             for i = 1:length(rampStarts)
%                 % flip logical for this ramp
%                 rampOnLogical(rampStarts(i):rampEnds(i)) = true;
%             end
% 
%             
%             % get ramp velocities, in degrees/sec
%             % volts per pixel
%             vPerPx = V_RANGE / pattern.x_num;
% 
%             % preallocate
%             rampVel = zeros(size(rampDurs));
% 
%             % loop through each ramp, get velocity
%             for i = 1:length(rampDurs)
%                 thisRampStart = rampStarts(i);
%                 thisRampEnd = rampEnds(i);
% 
%                 % max voltage and index of max voltage
%                 [maxV, maxVInd] = max(...
%                     daqData.panelsDAC0X(thisRampStart:thisRampEnd));
% 
%                 % number of LEDs moved
%                 numLEDs = round((maxV - spaceV) / vPerPx);
% 
%                 % speed of this ramp
%                 thisRampSpd = (numLEDs * DEG_PER_PX) / rampDurs(i);
% 
%                 % get direction of movement
%                 relStart = 1;
%                 relEnd = thisRampEnd - thisRampStart + 1;
% 
%                 % if max is closer to start than end of ramp
%                 %  stimulus rotating CW -> (+)
%                 if (abs(maxVInd - relStart) < abs(maxVInd - relEnd))
%                     rampVel(i) = thisRampSpd;
%                 % if max is closer to end than start   
%                 %  stimulus rotating CCW -> (-)
%                 else
%                     rampVel(i) = thisRampSpd * -1;
%                 end
%             end
% 
%             % convert actual ramp durations to commanded durations
%             cmdRampDurs = unique(inputParams.vInjParamsX.rampDur);
% 
%             rampsCmdDurs = zeros(size(rampDurs)); % preallocate
%         
%             for i = 1:length(rampDurs)
%                 thisRampDur = rampDurs(i);
%         
%                 % distance between this stimulation duration and all 
%                 %  commanded durations
%                 durDist = abs(cmdRampDurs - thisRampDur);
%                 % get the nearest commanded duration
%                 [~, nearInd] = min(durDist);
%                 thisCmdDur = cmdRampDurs(nearInd);
%         
%                 % update vector 
%                 rampsCmdDurs(i) = thisCmdDur;
%             end
% 
%             % convert actual ramp velocities to commanded velocities
%             cmdVel = (((inputParams.vInjParamsX.startV - ...
%                 inputParams.vInjParamsX.endV)/vPerPx) * DEG_PER_PX) ./ ...
%                 inputParams.vInjParamsX.rampDur;
% 
%             rampsCmdVels = zeros(size(rampDurs)); % preallocate
%             for i = 1:length(rampDurs)
%                 thisRampVel = rampVel(i);
%         
%                 % distance between this rampVel and all 
%                 %  commanded vel
%                 velDist = abs(cmdVel - thisRampVel);
%                 % get the nearest commanded velocity
%                 [~, nearInd] = min(velDist);
%                 thisCmdVel = cmdVel(nearInd);
%         
%                 % update vector 
%                 rampsCmdVels(i) = thisCmdVel;
%             end
% 
% 
%             % add to output struct
%             visstim.rampStartTimes = rampStartTimes;
%             visstim.rampEndTimes = rampEndTimes;
%             visstim.rampDurs = rampDurs;
%             visstim.rampCmdDurs = rampsCmdDurs;
%             visstim.rampOnLogical = rampOnLogical;
%             visstim.rampVels = rampVel;
%             visstim.rampCmdVels = rampsCmdVels;
%         end
% 
% 
%         % process Y channel signal - different depending on VInj protocol
%         % multiStepVInj
%         if (strcmp(inputParams.vInjProtocolY, 'multiStepVInj'))
%             % add visual stimulus Y protocol to output struct
%             visstim.visstimParams.vInjProtocolY = ...
%                 inputParams.vInjProtocolY;
% 
%             % number of possible voltage values
%             numSteps = pattern.y_num;
% 
%             % possible commanded output voltages
%             cmdOutV = (1:numSteps) * (V_RANGE/numSteps);
% 
%             % thresholds given possible output voltages
%             % set at midpoint between output voltages
%             stepVThresh = (cmdOutV(1:(end-1)) + cmdOutV(2:end))/2;
% 
%             % interpret steps given pattern
%             % regular square wave gratings
%             %  1 - grating, 2 - gray, 3 - off
%             if any(inputParams.patternIndex == REG_GRATING_PATTERN_IND)
%                 
%                 % get logical and start and end times for when grating
%                 %  displayed
%                 gratingLogical = (daqData.panelsDAC1Y' < stepVThresh(1));
% 
%                 % start and end indices for grating displayed
%                 gratingStarts = find(diff(gratingLogical) > 0.1);
%                 gratingEnds = find(diff(gratingLogical) < -0.1);
% 
%                 % start and end times for grating displayed
%                 gratingStartTimes = daqTime(gratingStarts);
%                 gratingEndTimes = daqTime(gratingEnds);
% 
%                 % get logical and start and end times for when gray
%                 %  displayed
%                 grayLogical = (daqData.panelsDAC1Y' > stepVThresh(1)) & ...
%                     (daqData.panelsDAC1Y' < stepVThresh(2));
% 
%                 % start and end indices for grating displayed
%                 grayStarts = find(diff(grayLogical) > 0.1);
%                 grayEnds = find(diff(grayLogical) < -0.1);
% 
%                 % start and end times for grating displayed
%                 grayStartTimes = daqTime(grayStarts);
%                 grayEndTimes = daqTime(grayEnds);
% 
%                 % add to output
%                 visstim.gratingStartTimes = gratingStartTimes;
%                 visstim.gratingEndTimes = gratingEndTimes;
%                 visstim.gratingOnLogical = gratingLogical;
%                 visstim.grayStartTimes = grayStartTimes;
%                 visstim.grayEndTimes = grayEndTimes;
%                 visstim.grayOnLogical = grayLogical;
%             end
%         end
% 
%         % for optomotor stimulus, combine ramp and grating times to get
%         %  static grating start and end
%         if any(inputParams.patternIndex == REG_GRATING_PATTERN_IND) && ...
%                 strcmp(inputParams.vInjProtocolX, 'spacedMultiRampVInj') && ...
%                 strcmp(inputParams.vInjProtocolY, 'multiStepVInj')
%             staticStartTimes = gratingStartTimes;
%             staticEndTimes = rampStartTimes;
%             staticLogical = gratingLogical & ~rampOnLogical;
% 
%             visstim.staticStartTimes = staticStartTimes;
%             visstim.staticEndTimes = staticEndTimes;
%             visstim.staticOnLogical = staticLogical;
%         end

%%
    % return empty struct if it does not meet criteria    
    else
        visstim = [];
    end

end