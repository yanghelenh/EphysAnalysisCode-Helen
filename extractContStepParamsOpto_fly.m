% extractContStepParamsOpto_fly.m
%
% Function to extract continuous leg step parameters with optogenetic 
%  stimulation.
% This uses the continuous estimate of the step parameters using the
%  envelope method (legStepsCont struct)
% Like extractLegStepParamsOpto_fly(), allows filtering valid opto trials
%  by the fly's behavior (has to be walking before, during, and
%  after the stimulation) beyond move/not move and selects no stim periods
%  as deliberate 'trials' during times of no opto (specifically, the 
%  middle).
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   durs - vector of all durations of stimulation to consider
%   NDs - vector of all NDs to consider
%   trialWindow - length 2 vector where 1st element is time before opto
%       starts to consider as trial and 2nd element is time after opto stim
%       ends to consider as part of the trial
%   walkTime - length 2 vector where 1st element is time before opto starts
%       and 2nd element is time after opto ends where the fly has to be
%       walking for the steps during the trial to be included
%   minWalkFwd - minimum average forward velocity fly needs to maintain
%       during window defined by walkTime for the trial to be included
%   flipLegsLR - boolean for whether to flip legs left right
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/6/23 - HHY
%
% UPDATED:
%   8/6/23 - HHY
%   8/18/23 - HHY - fix bug in inverting values of step parameters when
%       flipLegsLR flagged
%
function extractContStepParamsOpto_fly(durs, NDs, trialWindow, walkTime, ...
    minWalkFwd, flipLegsLR, pDataPath, saveFileDir)

    NUM_LEGS = 6;
    lrLegInd = [4 5 6 1 2 3]; % indices when flipping legs left/right

    % names of all step parameters to save
    stepParamNames = {'AEPX', 'PEPX', 'AEPY', 'PEPY', 'stepLengthX', ...
        'stepLengthY', 'stepLength', 'stepDirection'};

    % all the step parameters where values need to be * -1 when flipping
    %  legs left right
    flipStepParams = {'stepLengthY', 'AEPY', 'PEPY', 'stepDirection'};

    % circular step parameters
    circStepParams = {'stepDirection'};
    
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


    % preallocate
    for i = 1:length(stepParamNames)
        legStepsContOptoOne.(stepParamNames{i}) = [];
    end
    % to keep track of which ND for each trial
    legStepsContOptoOne.whichND = [];

    % make struct array, 1 for each duration
    legStepsContOpto = repmat(legStepsContOptoOne, numDurs,1);


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

        % check if this pData file has legStepsCont, opto, 
        %  fictracProc structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'opto')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % load data
        load(pDataFullPath, 'legStepsCont', 'opto', 'fictracProc');

        % get trial length, in frames
        trialLength = zeros(numDurs, 1);
        for j = 1:length(trialLength)
            ifi = median(diff(legStepsCont.t));
            trialLengthTime = durs(j) + trialWindow(1) + trialWindow(2);
            trialLength(j) = ceil(trialLengthTime/ifi);
        end

        % get trial times
        for j = 1:numDurs
            trialTimes{j} = (0:(trialLength(j)-1)) * ifi - trialWindow(1);
        end


        % loop through each opto stim trial, extract trial window, check if
        %  fly is walking
        for j = 1:length(opto.stimStartTimes)
            % start and end times of window where fly needs to be walking
            thisWalkStartTime = opto.stimStartTimes(j) - walkTime(1);
            thisWalkEndTime = opto.stimEndTimes(j) + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracProc.t>=thisWalkStartTime) & ...
                (fictracProc.t<=thisWalkEndTime);
            thisFwdVel = fictracProc.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                % start time of trial
                thisTrialStartTime = opto.stimStartTimes(j) - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    legStepsCont.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(opto.stimCmdDurs(j) == durs);
                thisTrialEndInd = thisTrialStartInd + ...
                    trialLength(thisDurInd) - 1;

                
                if (thisTrialStartInd <= length(legStepsCont.t)) && ...
                        (thisTrialEndInd <= length(legStepsCont.t))
                    % add this trial's step param to running tracker of trials
                    for k = 1:length(stepParamNames)
                        % flip leg index assignments
                         if (flipLegsLR)
                            % those parameters that need to be value inverted
                            if any(strcmpi(stepParamNames{k}, flipStepParams))
                                legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd) * -1);
                            else
                                legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd));
                            end
                         else
                            legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                legStepsCont.(stepParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd,:));
                         end
                    end
                    % add this trial's ND
                    legStepsContOpto(thisDurInd).whichND = [...
                        legStepsContOpto(thisDurInd).whichND; ...
                        opto.stimParams.ndFilter];
                end
            end
        end

        % loop through each no stim 'trial', defined by stim end time to
        %  stim start time
        for j = 1:(length(opto.stimStartTimes)-1)
            noStimDur = opto.stimStartTimes(j+1) - opto.stimEndTimes(j);
            thisStimDur = opto.stimCmdDurs(j);

            % walk time of this trial, as if there were actual stim
            % use stim duration of real stim of same index
            % 'stim' centered during no stim period
            thisWalkStartTime = opto.stimEndTimes(j) + noStimDur/2 - ...
                thisStimDur/2;
            % end time
            thisWalkEndTime = thisWalkStartTime + thisStimDur;

            % window where fly needs to be walking
            thisWalkStartTime = thisWalkStartTime - walkTime(1);
            thisWalkEndTime = thisWalkEndTime + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracProc.t>=thisWalkStartTime) & ...
                (fictracProc.t<=thisWalkEndTime);
            thisFwdVel = fictracProc.fwdVel(ftLog);

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                % start time of trial
                thisTrialStartTime = opto.stimEndTimes(j) + noStimDur/2 - ...
                    thisStimDur/2 - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    legStepsCont.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(opto.stimCmdDurs(j) == durs);
                thisTrialEndInd = thisTrialStartInd + ...
                    trialLength(thisDurInd) - 1;


                if (thisTrialStartInd <= length(legStepsCont.t)) && ...
                        (thisTrialEndInd <= length(legStepsCont.t))
                    % add this trial's step param to running tracker of trials
                    for k = 1:length(stepParamNames)
                        % flip leg index assignments
                         if (flipLegsLR)
                            % those parameters that need to be value inverted
                            if any(strcmpi(stepParamNames{k}, flipStepParams))
                                legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd) * -1);
                            else
                                legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd));
                            end
                         else
                            legStepsContOpto(thisDurInd).(stepParamNames{k}) = cat(3,...
                                legStepsContOpto(thisDurInd).(stepParamNames{k}), ...
                                legStepsCont.(stepParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd,:));
                         end
                    end
                    % add this trial's ND
                    legStepsContOpto(thisDurInd).whichND = [...
                        legStepsContOpto(thisDurInd).whichND; -1];
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % preallocate
    for i = 1:length(stepParamNames)
        for j = 1:numDurs
            for k = 1:numNDs
                legStepsContOptoMean(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
                legStepsContOptoSEM(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
                legStepsContOptoStdDev(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
            end
        end
    end
    
    % loop through all durations
    for i = 1:numDurs
        % loop through all NDs
        for j = 1:numNDs
            % get indices of trials that belong to this ND
            thisNDInd = find(legStepsContOpto(i).whichND == NDs(j));
            % check that there is data
            if ~isempty(thisNDInd)
                % loop through all step params
                for k = 1:length(stepParamNames)
                    % if circular
                    if any(strcmpi(stepParamNames{k}, circStepParams))
                        legStepsContOptoMean(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_mean(deg2rad(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd)),[], 3));
                        legStepsContOptoStdDev(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_std(deg2rad(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd)),[],[],3));
                        legStepsContOptoSEM(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_std(deg2rad(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd)),[], [],3)) / ...
                            sqrt(length(thisNDInd));
                    % if not circular    
                    else
                        legStepsContOptoMean(i,j).(stepParamNames{k}) = ...
                            mean(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd),3);
                        legStepsContOptoStdDev(i,j).(stepParamNames{k}) = ...
                            std(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd),[],3);
                        legStepsContOptoSEM(i,j).(stepParamNames{k}) = ...
                            std(legStepsContOpto(i).(stepParamNames{k})(:,:,thisNDInd),[],3) / ...
                            sqrt(length(thisNDInd));
                    end
                end
            end
        end
    end


    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_legStepsContOpto.mat'];
    save(saveFileFullName, 'legStepsContOpto', 'legStepsContOptoMean', ...
        'legStepsContOptoStdDev', 'legStepsContOptoSEM', 'trialWindow', ...
        'walkTime', 'minWalkFwd', 'flipLegsLR', 'durs', 'NDs', ...
        'trialTimes', '-v7.3');

end