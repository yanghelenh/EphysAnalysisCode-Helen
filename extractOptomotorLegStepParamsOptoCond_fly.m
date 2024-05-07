% extractOptomotorLegStepParamsOptoCond_fly.m
%
% Function to extract leg step parameters with optomotor visual stimuli
%  presented and optogenetic stimulation, with conditioning.
% This is a variant of extractLegStepParamsOptoCond_fly() with the addition
%  of optomotor visual stimuli. Works like conditioning
%  in extractOptomotorFicTracOptoCond_fly()
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   vels - vector of all optomotor stimulus velocities to consider. If 
%       vels = 0, combine all velocities
%   NDs - vector of all NDs to consider
%   trialWindow - length 2 vector where 1st element is time before
%       optomotor movement starts to consider as trial and 2nd element is 
%       time after optomotor movement stops to consider as part of the trial
%   cond - struct of conditioning, empty vector for no conditioning
%       whichParam - cell array (even if one parameter) of FicTrac
%           parameter names to condition on
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
%   flipLegsLR - boolean for whether to flip legs left right
%   pDataPath - path to folder containing pData files
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 1/28/24 - HHY
%
% UPDATED:
%   1/28/24 - HHY
%   4/12/24 - HHY - add conditioning on whether legs walking forward of
%       backwards
%   4/12/24 - HHY - remove outliers before computing mean, for non-circular
%       parameters
%   4/19/24 - HHY - conditioning changed to fictracSmo instead of
%       fictracProc, add pDataFNames
%   4/22/24 - HHY - add invAll to conditioning, change conditioning back to
%       fictracProc
%   4/29/24 - HHY - if vels = 0, trials are not split by vel and all are
%       assigned vel value of 0, regardless of actual vel
%
function extractOptomotorLegStepParamsOptoCond_fly(vels, NDs, trialWindow, ...
    cond, flipLegsLR, pDataPath, pDataFNames, saveFileDir)

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
    numConds = numNDs * numVels; 

    % key to map conditions to indices
    condKeyNDs = zeros(numConds, 1);
    condKeyVels = zeros(numVels, 1);
    counter = 1;
    for i = 1:numNDs
        for j = 1:numVels
            condKeyNDs(counter) = NDs(i);
            condKeyVels(counter) = vels(j);

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

    % add stepWhichLeg, visVel, optoND to All
    legStepsOptoAll.stance.stepWhichLeg = [];
    legStepsOptoAll.swing.stepWhichLeg = [];
    legStepsOptoAll.stance.visVel = [];
    legStepsOptoAll.swing.visVel = [];
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
        %  fictracSmo structs, if not, skip
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
        load(pDataFullPath, 'legSteps', 'opto', 'visstim', ...
            'fictracProc', 'fictracSmo');

        

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
        % velocity of valid trials
        trialVel = [];
        % ND of valid trials
        trialND = [];

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
%                     ftLog = (fictracProc.t>=thisCondStartTime) & ...
%                         (fictracProc.t<=thisCondEndTime);
                    ftLog = (fictracSmo.t>=thisCondStartTime) & ...
                        (fictracSmo.t<=thisCondEndTime);

                    % condition on step walking forward/backwards
                    if (strcmpi(cond.whichParam{k}, 'stepFwdBool'))
%                         thisTime = fictracProc.t(ftLog);
                        thisTime = fictracSmo.t(ftLog);
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
                        % condition on FicTracProc
%                         % the FicTrac parameter to condition on
%                         thisCondParam = fictracProc.(cond.whichParam{k});
%     
%                         % FicTrac parameter values during this time
%                         ftLog = (fictracProc.t>=thisCondStartTime) & ...
%                             (fictracProc.t<=thisCondEndTime);
                        
                        % condition on FicTracSmo
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

            % if this trial is valid (meetsCond true), record trial
            %  info
            if (meetsCond)
                % this trial start and end times, mod by trialWindow
                thisTrialStartTime = visstim.rampStartTimes(j) - ...
                    trialWindow(1);
                thisTrialEndTime = visstim.rampEndTimes(j) + ...
                    trialWindow(2);

                trialInd = [trialInd; j];
                trialStartTime = [trialStartTime; thisTrialStartTime];
                trialEndTime = [trialEndTime; thisTrialEndTime];
                % if vels = 0, then velocity is assigned value of 0,
                %  regardless of actual
                if (vels == 0)
                    trialVel = [trialVel; 0];
                else
                    trialVel = [trialVel; visstim.rampCmdVels(j)];
                end
                
                % this trial's ND
                % add this trial's ND or -1 if opto not on
                % check if opto was on during trial (during ramp)

                startOptoInd = find(...
                    visstim.rampStartTimes(j) >= opto.stimStartTimes,1,'last');
                endOptoInd = find(...
                    visstim.rampEndTimes(j) <= opto.stimEndTimes,1,'first');

                % with opto
                if (startOptoInd == endOptoInd)
                    trialND = [trialND; opto.stimParams.ndFilter];
                else % no opto
                    trialND = [trialND; -1];
                end
            end
        end

        % loop through all steps, add them to output vectors if they fall
        %  during valid trial times
        for j = 1:length(legSteps.stepWhichLeg)
            % loop through 2 half steps
            for k = 1:size(legSteps.stepLengths, 2)
                
                % this half step, start and end times
                startTime = legSteps.stepT(j,k);
                endTime = legSteps.stepT(j,k+1);

                % check if this step falls during valid trials
                % step has to start during stim window but can end after
                % start of valid windows 
                % index of trial this step falls into, if it does
                thisTrialInd = find((startTime >= trialStartTime) & ...
                    (startTime <= trialEndTime));
                % if this step falls during a trial, get step param and
                %  trial info and save into outputs
                if ~isempty(thisTrialInd)
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

                        % opto stim and vis vel params to add to step
                        legStepsOptoAll.stance.visVel = [...
                            legStepsOptoAll.stance.visVel; ...
                            trialVel(thisTrialInd)];
                        legStepsOptoAll.stance.optoND = [...
                            legStepsOptoAll.stance.optoND; ...
                            trialND(thisTrialInd)];

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

                        % opto stim and vis vel params to add to step
                        legStepsOptoAll.swing.visVel = [...
                            legStepsOptoAll.swing.visVel; ...
                            trialVel(thisTrialInd)];
                        legStepsOptoAll.swing.optoND = [...
                            legStepsOptoAll.swing.optoND; ...
                            trialND(thisTrialInd)];
                    end
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % loop through all conditions
    for i = 1:numConds
        thisND = condKeyNDs(i);
        thisVel = condKeyVels(i);

        % loop through all legs
        for j = 1:NUM_LEGS
            % loop through all step params
            for k = 1:length(stepParamNames)
                thisStanceValAll = legStepsOptoAll.stance.(stepParamNames{k});
                thisStanceVal = thisStanceValAll(...
                    (legStepsOptoAll.stance.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsOptoAll.stance.visVel == thisVel) & ...
                    (legStepsOptoAll.stance.optoND == thisND));

                thisSwingValAll = legStepsOptoAll.swing.(stepParamNames{k});
                thisSwingVal = thisSwingValAll(...
                    (legStepsOptoAll.swing.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsOptoAll.swing.visVel == thisVel) & ...
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
%                         % remove outliers
%                         thisStanceVal = rmoutliers(thisStanceVal);

                        legStepsOptoMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisStanceVal);
                        legStepsOptoStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal);
                        legStepsOptoSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)
%                         % remove outliers
%                         thisSwingVal = rmoutliers(thisSwingVal);

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
    saveFileFullName = [saveFileDir filesep flyName '_VisstimLegStepsOpto.mat'];
    save(saveFileFullName, 'legStepsOptoAll', 'legStepsOptoMeans', ...
        'legStepsOptoStdDev', 'legStepsOptoSEM', 'condKeyVels', ...
        'condKeyNDs', 'trialWindow', 'cond', ...
        'flipLegsLR', '-v7.3');

end