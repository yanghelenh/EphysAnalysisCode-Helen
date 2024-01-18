% extractLegStepParamsOptoCond_fly.m
%
% Function to extract leg step parameters with optogenetic stimulation.
% This is a variant of extractLegStepParamsOpto_fly() in which conditioning
%  on forward walking before/during/after the stimulation has been replaced
%  with a general user-specified set of conditions. Works like conditioning
%  in extractFicTracOptoCond_fly()
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   durs - vector of all durations of stimulation to consider
%   NDs - vector of all NDs to consider
%   optoTime - length 2 vector where 1st element is time after opto starts
%       to begin counting step as during opto (as time in sec relative to
%       opto start time) and 2nd element is time before opto ends to stop
%       counting step as during opto (as time in sec relative to opto end
%       time). If negative, start before opto stim starts/end after opto
%       stim ends.
%   cond - struct of conditioning, empty vector for no conditioning
%       whichParam - cell array (even if one parameter) of FicTrac
%           parameter names to condition on
%       cond - cell array (even if one condition) of conditions, matches up
%           to whichParam, as evaluatable strings
%       timeWin - cell array of vectors, matches up to whichParam and cond,
%           of time relative to opto stim where cond has to hold true.
%           Express each one as 4 element vector, 1 value for start and 1 
%           value for end, NaN for other 2: 
%           [start relative to opto start, start relative to opto end, 
%           end relative to opto start, end relative to opto end]
%           Negative for time before, positive for time after
%       during window defined by walkTime for the trial to be included
%   flipLegsLR - boolean for whether to flip legs left right
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 12/15/23 - HHY
%
% UPDATED:
%   12/17/23 - HHY
%
function extractLegStepParamsOptoCond_fly(durs, NDs, optoTime, cond, ...
    flipLegsLR, pDataPath, saveFileDir)

    NUM_LEGS = 6;

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY'};

    % all the step parameters where values need to be * -1 when flipping
    %  legs left right
    flipStepParams = {'stepYLengths', 'stepVelY', 'stepAEPY', ...
        'stepPEPY', 'stepDirections'};

    % circular step parameters
    circStepParams = {'stepDirections'};
    
    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
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
    % number of conditions (NDs by durs)
    numNDs = length(NDs);
    numDurs = length(durs);
    numConds = numNDs * numDurs; 

    % key to map conditions to indices
    condKeyNDs = zeros(numConds, 1);
    condKeyDurs = zeros(numDurs, 1);
    counter = 1;
    for i = 1:numNDs
        for j = 1:numDurs
            condKeyNDs(counter) = NDs(i);
            condKeyDurs(counter) = durs(j);

            counter = counter + 1;
        end
    end


    % preallocate
    for i = 1:length(stepParamNames)
        legStepsOptoAll.stance.(stepParamNames{i}) = [];

        legStepsOptoMeans.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsOptoStdDev.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsOptoSEM.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);

        legStepsOptoAll.swing.(stepParamNames{i}) = [];
        legStepsOptoMeans.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsOptoStdDev.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsOptoSEM.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
    end

    % add stepWhichLeg, optoDur, optoND to All
    legStepsOptoAll.stance.stepWhichLeg = [];
    legStepsOptoAll.swing.stepWhichLeg = [];
    legStepsOptoAll.stance.optoDur = [];
    legStepsOptoAll.swing.optoDur = [];
    legStepsOptoAll.stance.optoND = [];
    legStepsOptoAll.swing.optoND = [];

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

        % check if this pData file has legSteps, opto, 
        %  fictracProc structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'opto')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')))
            continue;
        end

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % load data
        load(pDataFullPath, 'legSteps', 'opto', 'fictracProc');

        

        % get matching b/w corresponding left and right legs
        rightLegInd = find(contains(legSteps.legIDs.names, 'R'));
        leftLegInd = find(contains(legSteps.legIDs.names, 'L'));
        matchedLegInd = zeros(length(rightLegInd),2);

        for j = 1:length(rightLegInd)
            thisLegNum = legSteps.legIDs.names{rightLegInd(j)}(end);
            thisLeftInd = find(contains(...
                legSteps.legIDs.names(leftLegInd),thisLegNum));
            matchedLegInd(j,1) = rightLegInd(j);
            matchedLegInd(j,2) = leftLegInd(thisLeftInd);
        end

        % allocate tracking for each trial
        % indices of valid trials
        trialInd = []; 
        % start time of valid trials, include mod by optoTime
        trialStartTime = [];
        % end time of valid trials, include mod by optoTime
        trialEndTime = [];
        % duration of valid trials
        trialDur = [];

        % loop through each opto stim trial, check if it meets walking
        %  criteria
        for j = 1:length(opto.stimStartTimes)
            % if there are no conditions, meets conditions by default
            if isempty(cond)
                % boolean for whether trial meets conditions
                meetsCond = true;
            else
                % initialize
                meetsCond = true;

                % loop through all conditions
                for k = 1:length(cond.whichParam)
                    % the FicTrac parameter to condition on
                    thisCondParam = fictracProc.(cond.whichParam{k});

                    % get start and end times of condition window
                    thisTimeWin = cond.timeWin{k};
                    % if first one is NaN, use second for start time
                    if isnan(thisTimeWin(1))
                        % second val for start time is relative to stim end
                        thisCondStartTime = opto.stimEndTimes(j) + ...
                            thisTimeWin(2);
                    else
                        % first val for start time is relative to stim
                        %  start
                        thisCondStartTime = opto.stimStartTimes(j) + ...
                            thisTimeWin(1);
                    end
                    % if third val is NaN, use fourth for end time
                    if isnan(thisTimeWin(3))
                        % fourth val for end time is relative to stim end
                        thisCondEndTime = opto.stimEndTimes(j) + ...
                            thisTimeWin(4);
                    else
                        % third val for end time is relative to stim start
                        thisCondEndTime = opto.stimStartTimes(j) + ...
                            thisTimeWin(3);
                    end

                    % FicTrac parameter values during this time
                    ftLog = (fictracProc.t>=thisCondStartTime) & ...
                        (fictracProc.t<=thisCondEndTime);
                    thisFTVal = thisCondParam(ftLog);

                    % check if criteria met
                    % logical for all time points
                    thisCondLog = eval(['thisFTVal' cond.cond{k}]);
                    % whether this condition is met (logical is all true)
                    thisCondMet = all(thisCondLog);

                    % combine with other conditions
                    % as soon as one condition is false, this will always
                    %  be false
                    meetsCond = meetsCond && thisCondMet; 
                end
            end

            % if this trial is valid (meetsCond true), record trial
            %  info
            if (meetsCond)
                trialInd = [trialInd; j];
                trialStartTime = [trialStartTime; ...
                    opto.stimStartTimes(j) + optoTime(1)];
                trialEndTime = [trialEndTime; ...
                    opto.stimEndTimes(j) - optoTime(2)];
                trialDur = [trialDur; opto.stimCmdDurs(j)];
            end
        end

        % loop through no stim periods, check if it meets walking criteria
        % allocate tracking for each no stim 'trial'
        % indices of valid trials
        nsTrialInd = []; 
        % start time of valid trials, include mod by optoTime
        nsTrialStartTime = [];
        % end time of valid trials, include mod by optoTime
        nsTrialEndTime = [];
        % duration of valid trials
        nsTrialDur = [];

        % loop through each no stim 'trial', defined by stim end time to
        %  stim start time
        for j = 1:(length(opto.stimStartTimes)-1)
            noStimDur = opto.stimStartTimes(j+1) - opto.stimEndTimes(j);
            thisStimDur = opto.stimCmdDurs(j);

            % this 'trial''s stim start and end time
            thisStimStartTime = opto.stimEndTimes(j) + noStimDur/2 - ...
                thisStimDur/2;
            % end time
            thisStimEndTime = thisStimStartTime + thisStimDur;

            % if there are no conditions, meets conditions by default
            if isempty(cond)
                % boolean for whether trial meets conditions
                meetsCond = true;
            else
                % initialize
                meetsCond = true;

                % loop through all conditions
                for k = 1:length(cond.whichParam)
                    % the FicTrac parameter to condition on
                    thisCondParam = fictracProc.(cond.whichParam{k});

                    % get start and end times of condition window
                    thisTimeWin = cond.timeWin{k};
                    % if first one is NaN, use second for start time
                    if isnan(thisTimeWin(1))
                        % second val for start time is relative to stim end
                        thisCondStartTime = thisStimEndTime + ...
                            thisTimeWin(2);
                    else
                        % first val for start time is relative to stim
                        %  start
                        thisCondStartTime = thisStimStartTime + ...
                            thisTimeWin(1);
                    end
                    % if third val is NaN, use fourth for end time
                    if isnan(thisTimeWin(3))
                        % fourth val for end time is relative to stim end
                        thisCondEndTime = thisStimEndTime + ...
                            thisTimeWin(4);
                    else
                        % third val for end time is relative to stim start
                        thisCondEndTime = thisStimStartTime + ...
                            thisTimeWin(3);
                    end

                    % FicTrac parameter values during this time
                    ftLog = (fictracProc.t>=thisCondStartTime) & ...
                        (fictracProc.t<=thisCondEndTime);
                    thisFTVal = thisCondParam(ftLog);

                    % check if criteria met
                    % logical for all time points
                    thisCondLog = eval(['thisFTVal' cond.cond{k}]);
                    % whether this condition is met (logical is all true)
                    thisCondMet = all(thisCondLog);

                    % combine with other conditions
                    % as soon as one condition is false, this will always
                    %  be false
                    meetsCond = meetsCond && thisCondMet; 
                end
            end

            % if this trial is valid (meetsCond true), record trial
            %  info
            if (meetsCond)
                nsTrialInd = [nsTrialInd; j];
                nsTrialStartTime = [nsTrialStartTime; ...
                    thisStimStartTime + optoTime(1)];
                nsTrialEndTime = [nsTrialEndTime; ...
                    thisStimEndTime - optoTime(2)];
                nsTrialDur = [nsTrialDur; opto.stimCmdDurs(j)];
            end
        end

        % loop through all steps, add them to output vectors if they fall
        %  during valid trial times (or no stim 'trial' times)
        for j = 1:length(legSteps.stepWhichLeg)
            % loop through 2 half steps
            for k = 1:size(legSteps.stepLengths, 2)
                
                % this half step, start and end times
                startTime = legSteps.stepT(j,k);
                endTime = legSteps.stepT(j,k+1);

                % check if this step falls during valid trials/no stim
                %  trials
                % step has to start during stim window but can end after
                % start of valid windows 
                % index of trial this step falls into, if it does
                thisTrialInd = find((startTime >= trialStartTime) & ...
                    (startTime <= trialEndTime));
                thisNsTrialInd = find((startTime >= nsTrialStartTime) & ...
                    (startTime <= nsTrialEndTime));
                % if this step falls during a trial, get step param and
                %  trial info and save into outputs
                if (~isempty(thisTrialInd) || ~isempty(thisNsTrialInd))
                    % get whether this is a swing or stance half step
                    thisWhichPhase = legSteps.stepSwingStance(j,k);

                    if (thisWhichPhase == 1) % stance   
                        curLegInd = legSteps.stepWhichLeg(j);

                        % if flip left and right legs
                        if (flipLegsLR)
                            % flip left and right legs
                            % r for row index, c for column index
                            [r, c] = ind2sub(size(matchedLegInd),...
                                find(matchedLegInd == curLegInd));
                            % matched ind are across columns, flip columns
                            if (c==1)
                                c = 2;
                            else
                                c = 1;
                            end
                            % new leg index, swapping left and right
                            legStepsOptoAll.stance.stepWhichLeg = ...
                                [legStepsOptoAll.stance.stepWhichLeg; ...
                                matchedLegInd(r,c)];
                        else
                            legStepsOptoAll.stance.stepWhichLeg = ...
                                [legStepsOptoAll.stance.stepWhichLeg; ...
                                curLegInd];
                        end

                        % loop through all step parameters
                        for l = 1:length(stepParamNames)
                            % if we need to flip legs left/right
                            if (flipLegsLR)
                                % those parameters that need to be adjusted
                                if any(strcmpi(stepParamNames{l}, flipStepParams))
                                    legStepsOptoAll.stance.(stepParamNames{l}) = ...
                                        [legStepsOptoAll.stance.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k) * -1];

                                % those parameters that don't need to be
                                %  adjusted
                                else
                                    legStepsOptoAll.stance.(stepParamNames{l}) = ...
                                        [legStepsOptoAll.stance.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k)];
                                end
                            else
                                legStepsOptoAll.stance.(stepParamNames{l}) = ...
                                    [legStepsOptoAll.stance.(stepParamNames{l}); ...
                                    legSteps.(stepParamNames{l})(j,k)];
                            end
                        end

                        % opto params to add to step
                        % ND = -1 if no stim trial
                        if (~isempty(thisTrialInd)) % stim trial
                            legStepsOptoAll.stance.optoDur = [...
                                legStepsOptoAll.stance.optoDur; ...
                                trialDur(thisTrialInd)];
                            legStepsOptoAll.stance.optoND = [...
                                legStepsOptoAll.stance.optoND; ...
                                opto.stimParams.ndFilter];
                        else % not stim trial
                           legStepsOptoAll.stance.optoDur = [...
                                legStepsOptoAll.stance.optoDur; ...
                                nsTrialDur(thisNsTrialInd)];
                            legStepsOptoAll.stance.optoND = [...
                                legStepsOptoAll.stance.optoND; -1];
                        end

                    elseif (thisWhichPhase == -1) % swing
                        curLegInd = legSteps.stepWhichLeg(j);

                        % if flip left and right legs
                        if (flipLegsLR)
                            % flip left and right legs
                            % r for row index, c for column index
                            [r, c] = ind2sub(size(matchedLegInd),...
                                find(matchedLegInd == curLegInd));
                            % matched ind are across columns, flip columns
                            if (c==1)
                                c = 2;
                            else
                                c = 1;
                            end
                            % new leg index, swapping left and right
                            legStepsOptoAll.swing.stepWhichLeg = ...
                                [legStepsOptoAll.swing.stepWhichLeg; ...
                                matchedLegInd(r,c)];
                        else
                            legStepsOptoAll.swing.stepWhichLeg = ...
                                [legStepsOptoAll.swing.stepWhichLeg; ...
                                curLegInd];
                        end

                        % loop through all step parameters
                        for l = 1:length(stepParamNames)
                            % if we need to flip legs left/right
                            if (flipLegsLR)
                                % those parameters that need to be adjusted
                                if any(strcmpi(stepParamNames{l}, flipStepParams))
                                    legStepsOptoAll.swing.(stepParamNames{l}) = ...
                                        [legStepsOptoAll.swing.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k) * -1];

                                % those parameters that don't need to be
                                %  adjusted
                                else
                                    legStepsOptoAll.swing.(stepParamNames{l}) = ...
                                        [legStepsOptoAll.swing.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k)];
                                end
                            else
                                legStepsOptoAll.swing.(stepParamNames{l}) = ...
                                    [legStepsOptoAll.swing.(stepParamNames{l}); ...
                                    legSteps.(stepParamNames{l})(j,k)];
                            end
                        end

                        % opto params to add to step
                        % ND = -1 if no stim trial
                        if (~isempty(thisTrialInd)) % stim trial
                            legStepsOptoAll.swing.optoDur = [...
                                legStepsOptoAll.swing.optoDur; ...
                                trialDur(thisTrialInd)];
                            legStepsOptoAll.swing.optoND = [...
                                legStepsOptoAll.swing.optoND; ...
                                opto.stimParams.ndFilter];
                        else % not stim trial
                           legStepsOptoAll.swing.optoDur = [...
                                legStepsOptoAll.swing.optoDur; ...
                                nsTrialDur(thisNsTrialInd)];
                            legStepsOptoAll.swing.optoND = [...
                                legStepsOptoAll.swing.optoND; -1];
                        end
                    end
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % loop through all conditions
    for i = 1:numConds
        thisND = condKeyNDs(i);
        thisDur = condKeyDurs(i);

        % loop through all legs
        for j = 1:NUM_LEGS
            % loop through all step params
            for k = 1:length(stepParamNames)
                thisStanceValAll = legStepsOptoAll.stance.(stepParamNames{k});
                thisStanceVal = thisStanceValAll(...
                    (legStepsOptoAll.stance.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsOptoAll.stance.optoDur == thisDur) & ...
                    (legStepsOptoAll.stance.optoND == thisND));

                thisSwingValAll = legStepsOptoAll.swing.(stepParamNames{k});
                thisSwingVal = thisSwingValAll(...
                    (legStepsOptoAll.swing.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsOptoAll.swing.optoDur == thisDur) & ...
                    (legStepsOptoAll.swing.optoND == thisND));

                % circular stats if circular
                if any(strcmpi(stepParamNames{k}, circStepParams))
                    % check that there's data, otherwise, will leave value
                    %  as NaN
                    if ~isempty(thisStanceVal)
                        % convert to radians and back
                        thisStanceVal = deg2rad(thisStanceVal);
                        legStepsOptoMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_mean(thisStanceVal));
                        legStepsOptoStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisStanceVal));
                        legStepsOptoSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisStanceVal)) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)
                        % convert to radians and back
                        thisSwingVal = deg2rad(thisSwingVal);
                        legStepsOptoMeans.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_mean(thisSwingVal));
                        legStepsOptoStdDev.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisSwingVal));
                        legStepsOptoSEM.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisSwingVal)) / sqrt(length(thisSwingVal));
                    end

                % otherwise, regular stats   
                else
                    % check that there's data, otherwise, will leave value
                    %  as NaN
                    if ~isempty(thisStanceVal)
                        legStepsOptoMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisStanceVal);
                        legStepsOptoStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal);
                        legStepsOptoSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)
                        legStepsOptoMeans.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisSwingVal);
                        legStepsOptoStdDev.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisSwingVal);
                        legStepsOptoSEM.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisSwingVal) / sqrt(length(thisSwingVal));
                    end
                end
            end
        end
    end

    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_legStepsOpto.mat'];
    save(saveFileFullName, 'legStepsOptoAll', 'legStepsOptoMeans', ...
        'legStepsOptoStdDev', 'legStepsOptoSEM', 'condKeyDurs', ...
        'condKeyNDs', 'optoTime', 'cond', ...
        'flipLegsLR', '-v7.3');

end