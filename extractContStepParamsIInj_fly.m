% extractContStepParamsIInj_fly.m
%
% Function to extract continuous leg step parameters with current injection 
%  stimulation.
% This uses the continuous estimate of the step parameters using the
%  envelope method (legStepsCont struct)
% Like extractLegStepParamsIInj_fly(), allows filtering valid IInj trials
%  by the fly's behavior (has to be walking before, during, and
%  after the stimulation) beyond move/not move and selects no stim periods
%  as deliberate 'trials' during times of no stim (specifically, the 
%  middle). Also, allows filtering based on spike rate during stimulation.
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   trialWindow - length 2 vector where 1st element is time before iInj
%       starts to consider as trial and 2nd element is time after iInj stim
%       ends to consider as part of the trial
%   walkTime - length 2 vector where 1st element is time before iInj starts
%       and 2nd element is time after iInj ends where the fly has to be
%       walking for the steps during the trial to be included
%   minWalkFwd - minimum average forward velocity fly needs to maintain
%       during window defined by walkTime for the trial to be included
%   spikeRateCond - cell array of spike rate conditions, same length as
%       numAmps * numDurs + 1, that mean spike rate during iInj has to
%       achieve during stimulation for trial to be included
%   flipLegsLR - boolean for whether to flip legs left right
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/21/23 - HHY
%
% UPDATED:
%   8/21/23 - HHY
%
function extractContStepParamsIInj_fly(amps, durs, trialWindow, walkTime, ...
    minWalkFwd, spikeRateCond, flipLegsLR, pDataPath, saveFileDir)

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
        legStepsContIInjOne.(stepParamNames{i}) = [];
    end
    % to keep track of which amp for each trial
    legStepsContIInjOne.whichAmp = [];

    % make struct array, 1 for each duration
    legStepsContIInj = repmat(legStepsContIInjOne, numDurs,1);


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

        % check if this pData file has legStepsCont, iInj, 
        %  fictracProc, ephysSpikes structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'iInj')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            continue;
        end

        % load data
        load(pDataFullPath, 'legStepsCont', 'iInj', 'fictracProc',...
            'ephysSpikes');

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


        % loop through each iInj stim trial, extract trial window, check if
        %  fly is walking, check if spike rate within range
        for j = 1:length(iInj.startTimes)
            % start and end times of window where fly needs to be walking
            thisWalkStartTime = iInj.startTimes(j) - walkTime(1);
            thisWalkEndTime = iInj.endTimes(j) + walkTime(2);

            % fwd velocity during this time
            ftLog = (fictracProc.t>=thisWalkStartTime) & ...
                (fictracProc.t<=thisWalkEndTime);
            thisFwdVel = fictracProc.fwdVel(ftLog);

            % start and end times of window where spike rate is to be
            %  checked
            thisStimStartTime = iInj.startTimes(j);
            thisStimEndTime = iInj.endTimes(j);

            % this duration and amplitude, to get appropriate eval cond
            thisDur = iInj.durs(j);
            thisAmp = iInj.amps(j);

            thisEvalCond = spikeRateCond{(thisDur==condKeyDurs) & ...
                (thisAmp==condKeyAmps)};

            % get spike rate during this stim (as number of
            %  spikes/duration)
            thisStartInd = find(ephysSpikes.t >= thisStimStartTime, 1, ...
                'first');
            thisEndInd = find(ephysSpikes.t <= thisStimEndTime, 1, 'last');

            numSpikes = sum((ephysSpikes.startInd >= thisStartInd) & ...
                (ephysSpikes.startInd <= thisEndInd));
            stimDur = ephysSpikes.t(thisEndInd) - ...
                ephysSpikes.t(thisStartInd);
            thisSpikeRate = numSpikes/stimDur;

            % if this trial is valid (fwd vel not below min and spike rate 
            %  meets criteria, record trial info
            if (~any(thisFwdVel < minWalkFwd) && ...
                    eval(['thisSpikeRate' thisEvalCond]))
                % start time of trial
                thisTrialStartTime = iInj.startTimes(j) - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    legStepsCont.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(iInj.durs(j) == durs);
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
                                legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd) * -1);
                            else
                                legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd));
                            end
                         else
                            legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                legStepsCont.(stepParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd,:));
                         end
                    end
                    % add this trial's amp
                    legStepsContIInj(thisDurInd).whichAmp = [...
                        legStepsContIInj(thisDurInd).whichAmp; ...
                        iInj.amps(j)];
                end
            end
        end

        % loop through each no stim 'trial', defined by stim end time to
        %  stim start time
        for j = 1:(length(iInj.startTimes)-1)
            noStimDur = iInj.startTimes(j+1) - iInj.endTimes(j);
            thisStimDur = iInj.durs(j);

            % walk time of this trial, as if there were actual stim
            % use stim duration of real stim of same index
            % 'stim' centered during no stim period
            thisWalkStartTime = iInj.endTimes(j) + noStimDur/2 - ...
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


            % start and end times of window where spike rate is to be
            %  checked
            thisStimStartTime = iInj.endTimes(j) + noStimDur/2 - ...
                thisStimDur/2;
            thisStimEndTime = thisWalkStartTime + thisStimDur;

            % this duration and amplitude, to get appropriate eval cond
            thisDur = iInj.durs(j);
            thisAmp = 0;

            thisEvalCond = spikeRateCond{(thisDur==condKeyDurs) & ...
                (thisAmp==condKeyAmps)};

            % get spike rate during this stim (as number of
            %  spikes/duration)
            thisStartInd = find(ephysSpikes.t >= thisStimStartTime, 1, ...
                'first');
            thisEndInd = find(ephysSpikes.t <= thisStimEndTime, 1, 'last');

            numSpikes = sum((ephysSpikes.startInd >= thisStartInd) & ...
                (ephysSpikes.startInd <= thisEndInd));
            stimDur = ephysSpikes.t(thisEndInd) - ...
                ephysSpikes.t(thisStartInd);
            thisSpikeRate = numSpikes/stimDur;

            % if this trial is valid (fwd vel not below min and spike rate 
            %  meets criteria, record trial info
            if (~any(thisFwdVel < minWalkFwd) && ...
                    eval(['thisSpikeRate' thisEvalCond]))
                % start time of trial
                thisTrialStartTime = iInj.endTimes(j) + noStimDur/2 - ...
                    thisStimDur/2 - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    legStepsCont.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(iInj.durs(j) == durs);
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
                                legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd) * -1);
                            else
                                legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                    legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                    legStepsCont.(stepParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd,lrLegInd));
                            end
                         else
                            legStepsContIInj(thisDurInd).(stepParamNames{k}) = cat(3,...
                                legStepsContIInj(thisDurInd).(stepParamNames{k}), ...
                                legStepsCont.(stepParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd,:));
                         end
                    end
                    % add this trial's amp (0 for no stim)
                    legStepsContIInj(thisDurInd).whichAmp = [...
                        legStepsContIInj(thisDurInd).whichAmp; 0];
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % preallocate
    for i = 1:length(stepParamNames)
        for j = 1:numDurs
            for k = 1:numAmps
                legStepsContIInjMean(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
                legStepsContIInjSEM(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
                legStepsContIInjStdDev(j,k).(stepParamNames{i}) = nan(...
                    trialLength(j), NUM_LEGS);
            end
        end
    end
    
    % loop through all durations
    for i = 1:numDurs
        % loop through all amps
        for j = 1:numAmps
            % get indices of trials that belong to this amp
            thisAmpInd = find(legStepsContIInj(i).whichAmp == amps(j));
            % check that there is data
            if ~isempty(thisAmpInd)
                % loop through all step params
                for k = 1:length(stepParamNames)
                    % if circular
                    if any(strcmpi(stepParamNames{k}, circStepParams))
                        legStepsContIInjMean(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_mean(deg2rad(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd)),[], 3));
                        legStepsContIInjStdDev(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_std(deg2rad(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd)),[],[],3));
                        legStepsContIInjSEM(i,j).(stepParamNames{k}) = ...
                            rad2deg(circ_std(deg2rad(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd)),[], [],3)) / ...
                            sqrt(length(thisAmpInd));
                    % if not circular    
                    else
                        legStepsContIInjMean(i,j).(stepParamNames{k}) = ...
                            mean(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd),3);
                        legStepsContIInjStdDev(i,j).(stepParamNames{k}) = ...
                            std(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd),[],3);
                        legStepsContIInjSEM(i,j).(stepParamNames{k}) = ...
                            std(legStepsContIInj(i).(stepParamNames{k})(:,:,thisAmpInd),[],3) / ...
                            sqrt(length(thisAmpInd));
                    end
                end
            end
        end
    end


    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_legStepsContIInj.mat'];
    save(saveFileFullName, 'legStepsContIInj', 'legStepsContIInjMean', ...
        'legStepsContIInjStdDev', 'legStepsContIInjSEM', 'trialWindow', ...
        'walkTime', 'minWalkFwd', 'flipLegsLR', 'durs', 'amps', ...
        'trialTimes', 'spikeRateCond', '-v7.3');

end