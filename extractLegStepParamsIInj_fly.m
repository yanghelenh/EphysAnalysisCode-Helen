% extractLegStepParamsIInj_fly.m
%
% Function to extract leg step parameters with current injection.
%  This is an updated variant of sortLegStepsByIInj() and
%  saveLegStepParamByCond_fly(), that does basically the same as those two
%  functions together, except that it allows for filtering valid I inj
%  trials by the fly's behavior (has to be walking before, during, and
%  after the stimulation) beyond move/not move and selects no stim periods
%  as deliberate 'trials' during times of no IInj (specifically, the 
%  middle). Also, operates on all step parameters instead of one at a time.
% Option to remove outliers before computing mean/std dev/SEM, only on
%  non-circular parameters
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   iInjTime - length 2 vector where 1st element is time after iInj starts
%       to begin counting step as during iInj (as time in sec relative to
%       iInj start time) and 2nd element is time before iInj ends to stop
%       counting step as during iInj (as time in sec relative to iInj end
%       time)
%   walkTime - length 2 vector where 1st element is time before iInj starts
%       and 2nd element is time after iInj ends where the fly has to be
%       walking for the steps during the trial to be included
%   minWalkFwd - minimum average forward velocity fly needs to maintain
%       during window defined by walkTime for the trial to be included
%   cond - struct of parameters about selecting specific legStep parameter
%     to condition on. [] for no conditioning
%       whichParam - string of legStep parameter to condition on
%       whichPhase - string for which phase of legStep parameter to
%           condition on ('stance' or 'swing')
%       cond - cell array of strings to condition on, for eval(), of size
%           numDurs * numAmps + 1; order: each dur for each amp
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
% CREATED: 8/18/23 - HHY
%
% UPDATED:
%   8/21/23 - HHY
%   8/24/23 - HHY - fix circular stats (forgot to convert to radians)
%   8/26/23 - HHY - fix bug in inverting step parameter values when
%       flipping left/right
%   5/7/24 - HHY - add abs val version of stepXLengths, stepYLengths; add
%       option for outlier removal; minWalkFwd on fictracSmo; add
%       pDataFNames as input
%   5/8/24 - HHY - minWalkFwd back on fictracProc
%
function extractLegStepParamsIInj_fly(amps, durs, iInjTime, walkTime, ...
    minWalkFwd, cond, flipLegsLR, outThresh, pDataPath, pDataFNames, ...
    saveFileDir)

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

    % add 0 for no stim condition to amps
    if isrow(amps)
        amps = [0 amps];
    else
        amps = [0; amps];
    end
    % number of conditions (amps by durs)
    numAmps = length(amps);
    numDurs = length(durs);
    numConds = numAmps * numDurs; 

    % key to map conditions to indices
    condKeyAmps = zeros(numConds, 1);
    condKeyDurs = zeros(numDurs, 1);
    counter = 1;
    for i = 1:numAmps
        for j = 1:numDurs
            condKeyAmps(counter) = amps(i);
            condKeyDurs(counter) = durs(j);

            counter = counter + 1;
        end
    end


    % preallocate
    for i = 1:length(stepParamNames)
        legStepsIInjAll.stance.(stepParamNames{i}) = [];

        legStepsIInjMeans.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsIInjStdDev.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsIInjSEM.stance.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);

        legStepsIInjAll.swing.(stepParamNames{i}) = [];
        legStepsIInjMeans.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsIInjStdDev.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
        legStepsIInjSEM.swing.(stepParamNames{i}) = nan(numConds, ...
            NUM_LEGS);
    end

    % add stepWhichLeg, duration, amplitude of stim to All
    legStepsIInjAll.stance.stepWhichLeg = [];
    legStepsIInjAll.swing.stepWhichLeg = [];
    legStepsIInjAll.stance.dur = [];
    legStepsIInjAll.swing.dur = [];
    legStepsIInjAll.stance.amp = [];
    legStepsIInjAll.swing.amp = [];

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

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % check if this pData file has legSteps, opto, 
        %  fictracProc structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'iInj')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')))
            continue;
        end

        % load data
        load(pDataFullPath, 'legSteps', 'iInj', 'fictracProc', ...
            'stanceStepParams', 'swingStepParams');

        % get stepXLengthsAbs and stepYLengthsAbs
        legSteps.stepXLengthsAbs = abs(legSteps.stepXLengths);
        stanceStepParams.stepXLengthsAbs = abs(stanceStepParams.stepXLengths);
        swingStepParams.stepXLengthsAbs = abs(swingStepParams.stepXLengths);

        legSteps.stepYLengthsAbs = abs(legSteps.stepYLengths);
        stanceStepParams.stepYLengthsAbs = abs(stanceStepParams.stepYLengths);
        swingStepParams.stepYLengthsAbs = abs(swingStepParams.stepYLengths); 

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
        % start time of valid trials, include mod by iInjTime
        trialStartTime = [];
        % end time of valid trials, include mod by iInjTime
        trialEndTime = [];
        % duration of valid trials
        trialDur = [];
        % amplitude of valid trials
        trialAmp = [];

        % loop through each I inj stim trial, check if it meets walking
        %  criteria
        for j = 1:length(iInj.startTimes)
            % start and end times of window where fly needs to be walking
            thisStartTime = iInj.startTimes(j) - walkTime(1);
            thisEndTime = iInj.endTimes(j) + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracProc.t>=thisStartTime) & ...
                (fictracProc.t<=thisEndTime);
            thisFwdVel = fictracProc.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                trialInd = [trialInd; j];
                trialStartTime = [trialStartTime; ...
                    iInj.startTimes(j) + iInjTime(1)];
                trialEndTime = [trialEndTime; ...
                    iInj.endTimes(j) - iInjTime(2)];
                trialDur = [trialDur; iInj.durs(j)];
                trialAmp = [trialAmp; iInj.amps(j)];
            end
        end

        % loop through no stim periods, check if it meets walking criteria
        % allocate tracking for each no stim 'trial'
        % indices of valid trials
        nsTrialInd = []; 
        % start time of valid trials, include mod by iInjTime
        nsTrialStartTime = [];
        % end time of valid trials, include mod by iInjTime
        nsTrialEndTime = [];
        % duration of valid trials
        nsTrialDur = [];
        % amplitude of valid trials
        nsTrialAmp = [];

        % loop through each no stim 'trial', defined by stim end time to
        %  stim start time
        for j = 1:(length(iInj.startTimes)-1)
            noStimDur = iInj.startTimes(j+1) - iInj.endTimes(j);
            thisStimDur = iInj.durs(j);

            % start time of this trial, as if there were actual stim
            % use stim duration of real stim of same index
            % 'stim' centered during no stim period
            thisTrialStartTime = iInj.endTimes(j) + noStimDur/2 - ...
                thisStimDur/2;
            % end time
            thisTrialEndTime = thisTrialStartTime + thisStimDur;

            % window where fly needs to be walking
            thisStartTime = thisTrialStartTime - walkTime(1);
            thisEndTime = thisTrialEndTime + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracProc.t>=thisStartTime) & ...
                (fictracProc.t<=thisEndTime);
            thisFwdVel = fictracProc.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                nsTrialInd = [nsTrialInd; j];
                nsTrialStartTime = [nsTrialStartTime; ...
                    thisTrialStartTime + iInjTime(1)];
                nsTrialEndTime = [nsTrialEndTime; ...
                    thisTrialEndTime - iInjTime(2)];
                nsTrialDur = [nsTrialDur; iInj.durs(j)];
                nsTrialAmp = [nsTrialAmp; 0];
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


                    % check if this step meets the cond criteria, when present
                    if (~isempty(cond))
                        % get amplitude and duration
                        % stim trial or no stim trial
                        if isempty(thisTrialInd) % no stim trial
                            thisAmp = nsTrialAmp(thisNsTrialInd);
                            thisDur = nsTrialDur(thisNsTrialInd);
                        else
                            thisAmp = trialAmp(thisTrialInd);
                            thisDur = trialDur(thisTrialInd);
                        end
                        % use amplitude and duration and keys to get index
                        %  into cond
                        thisCondLog = (condKeyAmps == thisAmp) & ...
                            (condKeyDurs == thisDur);
                        thisCond = cond.cond{thisCondLog};

                        % get value of conditioning parameter for this step
                        % swing or stance
                        if (strcmpi(cond.whichPhase, 'stance'))
                            thisCondVal = ...
                                stanceStepParams.(cond.whichParam)(j);
                        else
                            thisCondVal = ...
                                swingStepParams.(cond.whichParam)(j);
                        end
                        
                        % check if this step meets condition
                        % if not, go to next half step
                        if ~(eval(['thisCondVal' thisCond]))
                            continue;
                        end
                    end

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
                            legStepsIInjAll.stance.stepWhichLeg = ...
                                [legStepsIInjAll.stance.stepWhichLeg; ...
                                matchedLegInd(r,c)];
                        else
                            legStepsIInjAll.stance.stepWhichLeg = ...
                                [legStepsIInjAll.stance.stepWhichLeg; ...
                                curLegInd];
                        end

                        % loop through all step parameters
                        for l = 1:length(stepParamNames)
                            % if we need to flip legs left/right
                            if (flipLegsLR)
                                % those parameters that need to be adjusted
                                if any(strcmpi(stepParamNames{l}, flipStepParams))
                                    legStepsIInjAll.stance.(stepParamNames{l}) = ...
                                        [legStepsIInjAll.stance.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k) * -1];

                                % those parameters that don't need to be
                                %  adjusted
                                else
                                    legStepsIInjAll.stance.(stepParamNames{l}) = ...
                                        [legStepsIInjAll.stance.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k)];
                                end
                            else
                                legStepsIInjAll.stance.(stepParamNames{l}) = ...
                                    [legStepsIInjAll.stance.(stepParamNames{l}); ...
                                    legSteps.(stepParamNames{l})(j,k)];
                            end
                        end

                        % params to add to step
                        if (~isempty(thisTrialInd)) % stim trial
                            legStepsIInjAll.stance.dur = [...
                                legStepsIInjAll.stance.dur; ...
                                trialDur(thisTrialInd)];
                            legStepsIInjAll.stance.amp = [...
                                legStepsIInjAll.stance.amp; ...
                                trialAmp(thisTrialInd)];
                        else % not stim trial
                           legStepsIInjAll.stance.dur = [...
                                legStepsIInjAll.stance.dur; ...
                                nsTrialDur(thisNsTrialInd)];
                            legStepsIInjAll.stance.amp = [...
                                legStepsIInjAll.stance.amp; ...
                                nsTrialAmp(thisNsTrialInd)];
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
                            legStepsIInjAll.swing.stepWhichLeg = ...
                                [legStepsIInjAll.swing.stepWhichLeg; ...
                                matchedLegInd(r,c)];
                        else
                            legStepsIInjAll.swing.stepWhichLeg = ...
                                [legStepsIInjAll.swing.stepWhichLeg; ...
                                curLegInd];
                        end

                        % loop through all step parameters
                        for l = 1:length(stepParamNames)
                            % if we need to flip legs left/right
                            if (flipLegsLR)
                                % those parameters that need to be adjusted
                                if any(strcmpi(stepParamNames{l}, flipStepParams))
                                    legStepsIInjAll.swing.(stepParamNames{l}) = ...
                                        [legStepsIInjAll.swing.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k) * -1];

                                % those parameters that don't need to be
                                %  adjusted
                                else
                                    legStepsIInjAll.swing.(stepParamNames{l}) = ...
                                        [legStepsIInjAll.swing.(stepParamNames{l}); ...
                                        legSteps.(stepParamNames{l})(j,k)];
                                end
                            else
                                legStepsIInjAll.swing.(stepParamNames{l}) = ...
                                    [legStepsIInjAll.swing.(stepParamNames{l}); ...
                                    legSteps.(stepParamNames{l})(j,k)];
                            end
                        end

                        % params to add to step
                        if (~isempty(thisTrialInd)) % stim trial
                            legStepsIInjAll.swing.dur = [...
                                legStepsIInjAll.swing.dur; ...
                                trialDur(thisTrialInd)];
                            legStepsIInjAll.swing.amp = [...
                                legStepsIInjAll.swing.amp; ...
                                trialAmp(thisTrialInd)];
                        else % not stim trial
                           legStepsIInjAll.swing.dur = [...
                                legStepsIInjAll.swing.dur; ...
                                nsTrialDur(thisNsTrialInd)];
                            legStepsIInjAll.swing.amp = [...
                                legStepsIInjAll.swing.amp; ...
                                nsTrialAmp(thisNsTrialInd)];
                        end
                    end
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % loop through all conditions
    for i = 1:numConds
        thisAmp = condKeyAmps(i);
        thisDur = condKeyDurs(i);

        % loop through all legs
        for j = 1:NUM_LEGS
            % loop through all step params
            for k = 1:length(stepParamNames)
                thisStanceValAll = legStepsIInjAll.stance.(stepParamNames{k});
                thisStanceVal = thisStanceValAll(...
                    (legStepsIInjAll.stance.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsIInjAll.stance.dur == thisDur) & ...
                    (legStepsIInjAll.stance.amp == thisAmp));

                thisSwingValAll = legStepsIInjAll.swing.(stepParamNames{k});
                thisSwingVal = thisSwingValAll(...
                    (legStepsIInjAll.swing.stepWhichLeg == ...
                    legSteps.legIDs.ind(j)) & ...
                    (legStepsIInjAll.swing.dur == thisDur) & ...
                    (legStepsIInjAll.swing.amp == thisAmp));

                % circular stats if circular
                if any(strcmpi(stepParamNames{k}, circStepParams))
                    % check that there's data, otherwise, will leave value
                    %  as NaN
                    if ~isempty(thisStanceVal)
                        % convert to radians and back
                        thisStanceVal = deg2rad(thisStanceVal);
                        legStepsIInjMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_mean(thisStanceVal));
                        legStepsIInjStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisStanceVal));
                        legStepsIInjSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisStanceVal)) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)
                        % convert to radians and back
                        thisSwingVal = deg2rad(thisSwingVal);
                        legStepsIInjMeans.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_mean(thisSwingVal));
                        legStepsIInjStdDev.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            rad2deg(circ_std(thisSwingVal));
                        legStepsIInjSEM.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
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

                        legStepsIInjMeans.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisStanceVal);
                        legStepsIInjStdDev.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal);
                        legStepsIInjSEM.stance.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisStanceVal) / sqrt(length(thisStanceVal));
                    end

                    if ~isempty(thisSwingVal)

                        % outlier removal
                        if ~isempty(outThresh)
                            thisSwingVal = madRmOutliers(...
                                thisSwingVal, outThresh);
                        end

                        legStepsIInjMeans.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            mean(thisSwingVal);
                        legStepsIInjStdDev.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisSwingVal);
                        legStepsIInjSEM.swing.(stepParamNames{k})(i,legSteps.legIDs.ind(j)) = ...
                            std(thisSwingVal) / sqrt(length(thisSwingVal));
                    end
                end
            end
        end
    end

    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_legStepsIInj.mat'];
    save(saveFileFullName, 'legStepsIInjAll', 'legStepsIInjMeans', ...
        'legStepsIInjStdDev', 'legStepsIInjSEM', 'cond', 'condKeyDurs', ...
        'condKeyAmps', 'iInjTime', 'walkTime', 'minWalkFwd', ...
        'flipLegsLR', 'outThresh', '-v7.3');

end