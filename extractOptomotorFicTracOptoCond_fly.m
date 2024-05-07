% extractOptomotorFicTracOptoCond_fly.m
%
% Function to extract FicTrac values in response to optomotor visual 
%  stimulus with optogenetic stimulation.
% Adaptation of extractFicTracOptoCond_fly with addition of visual stimulus
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   vels - vector of all optomotor stimulus velocities to consider
%   NDs - vector of all optogenetic stim NDs to consider
%   trialWindow - length 2 vector where 1st element is time before
%       optomotor movement starts to consider as trial and 2nd element is 
%       time after optomotor movement stops to consider as part of the trial
%   cond - struct of conditioning, empty vector for no conditioning
%       whichParam - cell array (even if one parameter) of FicTrac
%           parameter names to condition on or 'stepFwdBool' for legs
%           stepping in forward walking direction
%       cond - cell array (even if one condition) of conditions, matches up
%           to whichParam, as evaluatable strings
%       legs - which leg (R1, R2, R3, L1, L2, L3) for stepFwdBool, [] for 
%           FicTrac
%       timeWin - cell array of vectors, matches up to whichParam and cond,
%           of time relative to optomotor visual stim where cond has to 
%           hold true.
%           Express each one as 4 element vector, 1 value for start and 1 
%           value for end, NaN for other 2: 
%           [start relative to opto start, start relative to opto end, 
%           end relative to opto start, end relative to opto end]
%           Negative for time before, positive for time after during window
%           defined by walkTime for the trial to be included
%       invAll - boolean for whether to condition on inversion of all
%           specified conditions ie: ~(A AND B) (not ~A AND ~B)
%   flipLR - boolean for whether to flip left/right asymmetric vars
%   pDataPath - path to folder containing pData files
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 1/22/24 - HHY
%
% UPDATED:
%   1/22/24 - HHY
%   1/28/24 - HHY - slight update to comments
%   1/29/24 - HHY - add conditioning on whether legs walking forward of
%       backwards
%   4/10/24 - HHY - add duration of optomotor stimulus, time of optogenetic
%       stimulation relative to optomotor as outputs
%   4/11/24 - HHY - fix bug where to round number of frames per trial
%       instead of using ceil, which was sensitive to very small changes in
%       ifi
%   4/22/24 - HHY - add invAll to conditioning, and pDataFNames to inputs
%   4/30/24 - HHY - if vels = 0, trials are not split by vel and all are
%       assigned vel value of 0, regardless of actual vel
%
function extractOptomotorFicTracOptoCond_fly(vels, NDs, trialWindow, ...
    cond, flipLR, pDataPath, pDataFNames, saveFileDir)

    % names of all fictrac parameters to save
    ftParamNames = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', ...
        'totAngSpd', 'totSpd', 'fwdCumPos', 'slideCumPos', ...
        'yawAngCumPos'};

    % all the FicTrac parameters where values need to be * -1 when flipping
    %  left right
    flipParams = {'slideVel', 'yawAngVel', 'slideCumPos', 'yawAngCumPos'};
    
    % prompt user to select pData files
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % add -1 for no stim condition to NDs
    if isrow(NDs)
        NDs = [-1 NDs];
    else
        NDs = [-1; NDs];
    end
    % number of conditions (NDs by vels)
    numNDs = length(NDs);
    numVels = length(vels);


    % preallocate
    for i = 1:length(ftParamNames)
        ftOptoOne.(ftParamNames{i}) = [];
    end
    % to keep track of which ND and which velocity for each trial
    ftOptoOne.whichND = [];
    ftOptoOne.whichVel = [];

    % make struct array, 1 for FicTrac values, 1 for change from start
    %  1 struct for each duration
    fictracOpto = repmat(ftOptoOne, numVels,1);
    changeFictracOpto = repmat(ftOptoOne, numVels,1);



    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has opto, visstim, fictracProc structs, 
        %  if not, skip
        if (~any(strcmpi(pDatVarsNames, 'opto')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'visstim')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % load data
        load(pDataFullPath, 'opto', 'fictracProc', 'visstim', ...
            'fictracSmo', 'legSteps');

        % get trial length, in frames
        % also, duration of optomotor during trial
        trialLength = zeros(numVels, 1);
        optomotorDur = zeros(numVels, 1);
        for j = 1:length(trialLength)
            ifi = median(diff(fictracProc.t));

            % if combining all velocities (vels = 0), assumes duration is
            %  same across all of them
            if (vels == 0)
                thisVelDur = median(visstim.rampCmdDurs);
            else
                thisVelDur = median(...
                    visstim.rampCmdDurs(visstim.rampCmdVels == vels(j)));
            end

            trialLengthTime = thisVelDur + trialWindow(1) + trialWindow(2);
            trialLength(j) = round(trialLengthTime/ifi);
            optomotorDur(j) = thisVelDur;
        end

        % pretrial window, in frames
        preTrialLength = floor(trialWindow(1)/ifi); 

        % get trial times
        for j = 1:numVels
            trialTimes{j} = (0:(trialLength(j)-1)) * ifi - trialWindow(1);
        end


        % loop through each optomotor  trial, check if it meets all cond
        %  criteria
        for j = 1:length(visstim.rampStartTimes)
            % if there are no conditions, meets conditions by default
            if isempty(cond)
                % boolean for whether trial meets conditions
                meetsCond = true;
            else
                % initialize
                meetsCond = true;

                % loop through all conditions
                for k = 1:length(cond.whichParam)
                    % get start and end times of condition window
                    thisTimeWin = cond.timeWin{k};
                    % if first one is NaN, use second for start time
                    if isnan(thisTimeWin(1))
                        % second val for start time is relative to stim end
                        thisCondStartTime = visstim.rampEndTimes(j) + ...
                            thisTimeWin(2);
                    else
                        % first val for start time is relative to stim
                        %  start
                        thisCondStartTime = visstim.rampStartTimes(j) + ...
                            thisTimeWin(1);
                    end
                    % if third val is NaN, use fourth for end time
                    if isnan(thisTimeWin(3))
                        % fourth val for end time is relative to stim end
                        thisCondEndTime = visstim.rampEndTimes(j) + ...
                            thisTimeWin(4);
                    else
                        % third val for end time is relative to stim start
                        thisCondEndTime = visstim.rampStartTimes(j) + ...
                            thisTimeWin(3);
                    end

                    % logical into fictrac for cond time
                    ftLog = (fictracProc.t>=thisCondStartTime) & ...
                        (fictracProc.t<=thisCondEndTime);

                    % condition on step walking forward/backwards
                    if (strcmpi(cond.whichParam{k}, 'stepFwdBool'))
                        thisTime = fictracProc.t(ftLog);
                        thisFwdLog = getStepFwdLogical(legSteps, thisTime,...
                            cond.legs{k});
                        if ~(eval(cond.cond{k})) % if target is false, invert
                            thisLog = ~thisFwdLog;
                        else
                            thisLog = thisFwdLog;
                        end
                        % if forward walking met for all time points
                        thisCondMet = all(thisLog);
                    else % condition on FicTrac parameter
%                         % condition on fictracProc
%                         % the FicTrac parameter to condition on
%                         thisCondParam = fictracProc.(cond.whichParam{k});
%     
%                         % FicTrac parameter values during this time
%                         ftLog = (fictracProc.t>=thisCondStartTime) & ...
%                             (fictracProc.t<=thisCondEndTime);

                        % condition on fictracSmo
                        % the FicTrac parameter to condition on
                        thisCondParam = fictracSmo.(cond.whichParam{k});
    
                        % FicTrac parameter values during this time
                        ftLog = (fictracSmo.t>=thisCondStartTime) & ...
                            (fictracSmo.t<=thisCondEndTime);

                        thisFTVal = thisCondParam(ftLog);
    
                        % check if criteria met
                        % logical for all time points
                        thisCondLog = eval(['thisFTVal' cond.cond{k}]);
                        % whether this condition is met (logical is all true)
                        thisCondMet = all(thisCondLog);
                    end
    
                    % combine with other conditions
                    % as soon as one condition is false, this will always
                    %  be false
                    meetsCond = meetsCond && thisCondMet; 
                end

                % if invert all conditioning
                if (cond.invAll)
                    meetsCond = ~meetsCond;
                end
            end

            % if this trial is valid (meets all cond criteria), record 
            %  trial info
            if (meetsCond)
                % start time of trial
                thisTrialStartTime = visstim.rampStartTimes(j) - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    fictracProc.t, 1, 'last');
                thisRampStartInd = find(visstim.rampStartTimes(j) >= ...
                    fictracProc.t, 1, 'last') ;
                thisRampEndInd = find(visstim.rampEndTimes(j) <= ...
                    fictracProc.t, 1, 'last') ;

                % get end index, based on velocity
                % if vels = 0, pool across velocities
                if (vels == 0)
                    thisVelInd = 1;
                else
                    thisVelInd = find(visstim.rampCmdVels(j) == vels);
                end

                % if this velocity is not to be considered, skip this trial
                if isempty(thisVelInd)
                    continue;
                end

                thisTrialEndInd = thisTrialStartInd + ...
                    trialLength(thisVelInd) - 1;
                % end index right before optomotor starts
                thisTrialPreEndInd = thisTrialStartInd + preTrialLength;


                % check whether any times during trial have dropped FicTrac
                thisTrialAllInd = thisTrialStartInd:thisTrialEndInd;
                overlapInd = intersect(thisTrialAllInd, fictracProc.dropInd);
                % if there is overlap, skip this trial
                if ~isempty(overlapInd)
                    continue;
                end

                
                if (thisTrialStartInd <= length(fictracProc.t)) && ...
                        (thisTrialEndInd <= length(fictracProc.t))
                    % add this trial's FicTrac to running tracker of trials
                    for k = 1:length(ftParamNames)
                        % flip parameter values if needed
                         if (flipLR)
                            % those parameters that need to be value inverted
                            if any(strcmpi(ftParamNames{k}, flipParams))
                                % FicTrac values
                                fictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisVelInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1);

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd) * -1);

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisVelInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1) - ...
                                    preOptoMean);

                            else
                                % FicTrac values
                                fictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisVelInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisVelInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preOptoMean);
                            end
                         else
                            fictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                fictracOpto(thisVelInd).(ftParamNames{k}), ...
                                fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd));
                            
                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));
                            
                                % mean subtracted FicTrac value
                                changeFictracOpto(thisVelInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisVelInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preOptoMean);
                         end
                    end
                    % add this trial's ND or -1 if opto not on
                    % check if opto was on during trial (during ramp)
                    thisRampStartTime = visstim.rampStartTimes(j);
                    thisRampEndTime = visstim.rampEndTimes(j);

                    startOptoInd = find(...
                        thisRampStartTime >= opto.stimStartTimes,1,'last');
                    endOptoInd = find(...
                        thisRampEndTime <= opto.stimEndTimes,1,'first');

                    if (startOptoInd == endOptoInd)
                        fictracOpto(thisVelInd).whichND = [...
                            fictracOpto(thisVelInd).whichND; ...
                            opto.stimParams.ndFilter];
                        changeFictracOpto(thisVelInd).whichND = [...
                            changeFictracOpto(thisVelInd).whichND; ...
                            opto.stimParams.ndFilter];
                    else % no opto
                        fictracOpto(thisVelInd).whichND = [...
                            fictracOpto(thisVelInd).whichND; ...
                            -1];
                        changeFictracOpto(thisVelInd).whichND = [...
                            changeFictracOpto(thisVelInd).whichND; ...
                            -1];
                    end

                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % preallocate
    for i = 1:length(ftParamNames)
        for j = 1:numVels
            for k = 1:numNDs
                % FicTrac values
                fictracOptoMean(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                fictracOptoSEM(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                fictracOptoStdDev(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);

                % change in FicTrac values
                changeFictracOptoMean(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                changeFictracOptoSEM(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                changeFictracOptoStdDev(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
            end
        end
    end
    
    % loop through all velocities
    for i = 1:numVels
        % loop through all NDs
        for j = 1:numNDs
            % get indices of trials that belong to this ND
            thisNDInd = find(fictracOpto(i).whichND == NDs(j));
            % check that there is data
            if ~isempty(thisNDInd)
                % loop through all FicTrac params
                for k = 1:length(ftParamNames)
                    % FicTrac values
                    fictracOptoMean(i,j).(ftParamNames{k}) = ...
                        mean(fictracOpto(i).(ftParamNames{k})(:,thisNDInd),2);
                    fictracOptoStdDev(i,j).(ftParamNames{k}) = ...
                        std(fictracOpto(i).(ftParamNames{k})(:,thisNDInd),[],2);
                    fictracOptoSEM(i,j).(ftParamNames{k}) = ...
                        std(fictracOpto(i).(ftParamNames{k})(:,thisNDInd),[],2) / ...
                        sqrt(length(thisNDInd));

                    % change in FicTrac values
                    changeFictracOptoMean(i,j).(ftParamNames{k}) = ...
                        mean(changeFictracOpto(i).(ftParamNames{k})(:,thisNDInd),2);
                    changeFictracOptoStdDev(i,j).(ftParamNames{k}) = ...
                        std(changeFictracOpto(i).(ftParamNames{k})(:,thisNDInd),[],2);
                    changeFictracOptoSEM(i,j).(ftParamNames{k}) = ...
                        std(changeFictracOpto(i).(ftParamNames{k})(:,thisNDInd),[],2) / ...
                        sqrt(length(thisNDInd));
                end
            end
        end
    end

    % get opto stim duration and start/end time (relative to optomotor
    %  stimulus)
    optoDur = opto.stimParams.allStimDurs;
    % get time before optomotor starts
    optomotorInd = find(visstim.visstimParams.yCmdAmps == 5);
    timeBfOptomotor = sum(visstim.visstimParams.yCmdDurs(1:(optomotorInd-1)));

    % get opto start
    optoStartTime = 0 - timeBfOptomotor + opto.stimParams.durBfStims;
    % get opto end
    optoEndTime = optoStartTime + optoDur;

    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_VisstimFicTracOpto.mat'];
    save(saveFileFullName, 'fictracOpto', 'fictracOptoMean', ...
        'fictracOptoStdDev', 'fictracOptoSEM', 'changeFictracOpto', ...
        'changeFictracOptoMean', 'changeFictracOptoStdDev', ...
        'changeFictracOptoSEM', 'trialWindow', 'optomotorDur', ...
        'optoDur', 'optoStartTime', 'optoEndTime', ...
        'cond', 'flipLR', 'vels', 'NDs', ...
        'trialTimes', '-v7.3');

end