% extractLegStepParamsOpto_fly.m
%
% Function to extract leg step parameters with optogenetic stimulation.
%  This is an updated variant of sortLegStepsByOpto() and
%  saveLegStepParamByCond_fly(), that does basically the same as those two
%  functions together, except that it allows for filtering valid opto
%  trials by the fly's behavior (has to be walking before, during, and
%  after the stimulation) beyond move/not move and selects no stim periods
%  as deliberate 'trials' during times of no opto (specifically, the 
%  middle). Also, operates on all step parameters instead of one at a time.
% Option to remove outliers before computing mean/std dev/SEM, only on
%  non-circular parameters
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
%   walkTime - length 2 vector where 1st element is time before opto starts
%       and 2nd element is time after opto ends where the fly has to be
%       walking for the steps during the trial to be included
%   minWalkFwd - minimum average forward velocity fly needs to maintain
%       during window defined by walkTime for the trial to be included
%   flipLegsLR - boolean for whether to flip legs left right
%   outThresh - threshold in number of MAD away to consider as outlier 
%       [] for no outlier removal 
%   pDataPath - path to folder containing pData files
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/4/23 - HHY
%
% UPDATED:
%   8/5/23 - HHY
%   8/21/23 - HHY - fix bug: stepYLengths also needs to be inverted
%       left/right
%   8/24/23 - HHY - fix circular stats: forgot to convert to radians
%   8/26/23 - HHY - fix bug in inverting step parameter values when
%       flipping left/right
%   5/7/24 - HHY - add abs val version of stepXLengths, stepYLengths; add
%       option for outlier removal; minWalkFwd on fictracSmo
%
function extractLegStepParamsOpto_fly(durs, NDs, optoTime, walkTime, ...
    minWalkFwd, flipLegsLR, outThresh, pDataPath, pDataFNames, saveFileDir)

    NUM_LEGS = 6;

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY', 'stepXLengthsAbs', 'stepYLengthsAbs'};

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
                ~any(strcmpi(pDatVarsNames, 'fictracSmo')))
            continue;
        end

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % load data
        load(pDataFullPath, 'legSteps', 'opto', 'fictracSmo');

        % get stepXLengthsAbs and stepYLengthsAbs
        legSteps.stepXLengthsAbs = abs(legSteps.stepXLengths);
        legSteps.stepYLengthsAbs = abs(legSteps.stepYLengths);

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
        
%         optoStartTime = [];

        % loop through each opto stim trial, check if it meets walking
        %  criteria
        for j = 1:length(opto.stimStartTimes)
            % start and end times of window where fly needs to be walking
            thisStartTime = opto.stimStartTimes(j) - walkTime(1);
            thisEndTime = opto.stimEndTimes(j) + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracSmo.t>=thisStartTime) & ...
                (fictracSmo.t<=thisEndTime);
            thisFwdVel = fictracSmo.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                trialInd = [trialInd; j];
                trialStartTime = [trialStartTime; ...
                    opto.stimStartTimes(j) + optoTime(1)];
%                 optoStartTime = [optoStartTime; opto.stimStartTimes(j)];
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

%         nsOptoStartTime = [];

        % loop through each no stim 'trial', defined by stim end time to
        %  stim start time
        for j = 1:(length(opto.stimStartTimes)-1)
            noStimDur = opto.stimStartTimes(j+1) - opto.stimEndTimes(j);
            thisStimDur = opto.stimCmdDurs(j);

            % start time of this trial, as if there were actual stim
            % use stim duration of real stim of same index
            % 'stim' centered during no stim period
            thisTrialStartTime = opto.stimEndTimes(j) + noStimDur/2 - ...
                thisStimDur/2;
            % end time
            thisTrialEndTime = thisTrialStartTime + thisStimDur;

            % window where fly needs to be walking
            thisStartTime = thisTrialStartTime - walkTime(1);
            thisEndTime = thisTrialEndTime + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracSmo.t>=thisStartTime) & ...
                (fictracSmo.t<=thisEndTime);
            thisFwdVel = fictracSmo.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                nsTrialInd = [nsTrialInd; j];
                nsTrialStartTime = [nsTrialStartTime; ...
                    thisTrialStartTime + optoTime(1)];
                nsTrialEndTime = [nsTrialEndTime; ...
                    thisTrialEndTime - optoTime(2)];

%                 nsOptoStartTime = [nsOptoStartTime; ...
%                     thisTrialStartTime];
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
%                 thisTrialInd = find((startTime >= trialStartTime) & ...
%                     (startTime < optoStartTime) & ...
%                     (endTime > optoStartTime) & ...
%                     (endTime <= trialEndTime));
%                 thisNsTrialInd = find((startTime >= nsTrialStartTime) & ...
%                     (startTime < nsOptoStartTime) & ...
%                     (endTime > nsOptoStartTime) & ...
%                     (endTime <= nsTrialEndTime));
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

                        % outlier removal
                        if ~isempty(outThresh)
                            thisStanceVal = madRmOutliers(...
                                thisStanceVal, outThresh);
                        end

                        legStepsOptoMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisStanceVal);
                        legStepsOptoStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal);
                        legStepsOptoSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)

                        % outlier removal
                        if ~isempty(outThresh)
                            thisSwingVal = madRmOutliers(...
                                thisSwingVal, outThresh);
                        end

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
        'condKeyNDs', 'optoTime', 'walkTime', 'minWalkFwd', ...
        'flipLegsLR', 'outThresh', '-v7.3');

end