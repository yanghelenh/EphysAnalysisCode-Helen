% extractFicTracIInj_fly.m
%
% Function to extract FicTrac values with current injection.
% Like extractContStepParamsIInj_fly(), allows filtering valid IInj trials
%  by the fly's behavior (has to be walking before, during, and
%  after the stimulation) beyond move/not move and selects no stim periods
%  as deliberate 'trials' during times of no IInj (specifically, the 
%  middle).
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   trialWindow - length 2 vector where 1st element is time before opto
%       starts to consider as trial and 2nd element is time after opto stim
%       eamps to consider as part of the trial
%   walkTime - length 2 vector where 1st element is time before opto starts
%       and 2nd element is time after opto eamps where the fly has to be
%       walking for the steps during the trial to be included
%   minWalkFwd - minimum average forward velocity fly needs to maintain
%       during window defined by walkTime for the trial to be included
%   spikeRateCond - cell array of spike rate conditions, same length as
%       numAmps * numDurs + 1, that mean spike rate during iInj has to
%       achieve during stimulation for trial to be included
%   flipLR - boolean for whether to flip left/right asymmetric vars
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/22/23 - HHY
%
% UPDATED:
%   8/22/23 - HHY
%
function extractFicTracIInj_fly(amps, durs, trialWindow, walkTime, ...
    minWalkFwd, spikeRateCond, flipLR, pDataPath, saveFileDir)

    % names of all fictrac parameters to save
    ftParamNames = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', ...
        'totAngSpd', 'totSpd', 'fwdCumPos', 'slideCumPos', ...
        'yawAngCumPos'};

    % all the FicTrac parameters where values need to be * -1 when flipping
    %  left right
    flipParams = {'slideVel', 'yawAngVel', 'slideCumPos', 'yawAngCumPos'};
    
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
    for i = 1:length(ftParamNames)
        ftIInjOne.(ftParamNames{i}) = [];
    end
    % to keep track of which amp for each trial
    ftIInjOne.whichAmp = [];

    % make struct array, 1 for FicTrac values, 1 for change from start
    %  1 struct for each duration
    fictracIInj = repmat(ftIInjOne, numDurs,1);
    changeFictracIInj = repmat(ftIInjOne, numDurs,1);



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

        % check if this pData file has opto, fictracProc structs, if not, 
        %  skip
        if (~any(strcmpi(pDatVarsNames, 'iInj')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            continue;
        end

        % load data
        load(pDataFullPath, 'iInj', 'fictracProc', 'ephysSpikes');

        % get trial length, in frames
        trialLength = zeros(numDurs, 1);
        for j = 1:length(trialLength)
            ifi = median(diff(fictracProc.t));
            trialLengthTime = durs(j) + trialWindow(1) + trialWindow(2);
            trialLength(j) = ceil(trialLengthTime/ifi);
        end

        % pretrial window, in frames
        preTrialLength = floor(trialWindow(1)/ifi); 

        % get trial times
        for j = 1:numDurs
            trialTimes{j} = (0:(trialLength(j)-1)) * ifi - trialWindow(1);
        end


        % loop through each IInj trial, extract trial window, check if
        %  fly is walking
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

            % if this trial is valid (fwd vel not below min), and spike rate 
            %  meets criteria, record trial info
            if (~any(thisFwdVel < minWalkFwd) && ...
                    eval(['thisSpikeRate' thisEvalCond]))
                % start time of trial
                thisTrialStartTime = iInj.startTimes(j) - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    fictracProc.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(iInj.durs(j) == durs);
                % if this duration is not to be considered, skip this trial
                if isempty(thisDurInd)
                    continue;
                end

                thisTrialEndInd = thisTrialStartInd + ...
                    trialLength(thisDurInd) - 1;
                % end index right before iInj starts
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
                                fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1);

                                % change in FicTrac values
                                % get mean before iInj
                                preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd) * -1);

                                % mean subtracted FicTrac value
                                changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1) - ...
                                    preIInjMean);

                            else
                                % FicTrac values
                                fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before iInj
                                preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preIInjMean);
                            end
                         else
                            fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before iInj
                                preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preIInjMean);
                         end
                    end
                    % add this trial's amp
                    fictracIInj(thisDurInd).whichAmp = [...
                        fictracIInj(thisDurInd).whichAmp; ...
                        iInj.amps(j)];
                    changeFictracIInj(thisDurInd).whichAmp = [...
                        changeFictracIInj(thisDurInd).whichAmp; ...
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

            % if this trial is valid (fwd vel not below min), record trial
            %  info
            if (~any(thisFwdVel < minWalkFwd))
                % start time of trial
                thisTrialStartTime = iInj.endTimes(j) + noStimDur/2 - ...
                    thisStimDur/2 - trialWindow(1);
                thisTrialStartInd = find(thisTrialStartTime >= ...
                    fictracProc.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(iInj.durs(j) == durs);
                
                % if this duration is not to be considered, skip this trial
                if isempty(thisDurInd)
                    continue;
                end

                thisTrialEndInd = thisTrialStartInd + ...
                    trialLength(thisDurInd) - 1;
                % end index right before opto starts
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
                    % add this trial's step param to running tracker of trials
                    for k = 1:length(ftParamNames)
                        % flip left/right
                         if (flipLR)
                            % those parameters that need to be value inverted
                            if any(strcmpi(ftParamNames{k}, flipParams))
                                % FicTrac values
                                fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1);

                                % change in FicTrac values
                                % get mean before iInj
                                preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd) * -1);

                                % mean subtracted FicTrac value
                                changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1) - ...
                                    preIInjMean);
                            else
                                % FicTrac values
                                fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before iInj
                                preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preIInjMean);
                            end
                         else
                            % FicTrac values
                            fictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                fictracIInj(thisDurInd).(ftParamNames{k}), ...
                                fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd));

                            % change in FicTrac values
                            % get mean before iInj
                            preIInjMean = mean(fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialPreEndInd));

                            % mean subtracted FicTrac value
                            changeFictracIInj(thisDurInd).(ftParamNames{k}) = cat(2,...
                                changeFictracIInj(thisDurInd).(ftParamNames{k}), ...
                                (fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd)) - ...
                                preIInjMean);
                         end
                    end
                    % add this trial's amp: 0 for no stim
                    fictracIInj(thisDurInd).whichAmp = [...
                        fictracIInj(thisDurInd).whichAmp; 0];
                    changeFictracIInj(thisDurInd).whichAmp = [...
                        changeFictracIInj(thisDurInd).whichAmp; 0];
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % preallocate
    for i = 1:length(ftParamNames)
        for j = 1:numDurs
            for k = 1:numAmps
                % FicTrac values
                fictracIInjMean(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                fictracIInjSEM(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                fictracIInjStdDev(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);

                % change in FicTrac values
                changeFictracIInjMean(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                changeFictracIInjSEM(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
                changeFictracIInjStdDev(j,k).(ftParamNames{i}) = nan(...
                    trialLength(j), 1);
            end
        end
    end
    
    % loop through all durations
    for i = 1:numDurs
        % loop through all amps
        for j = 1:numAmps
            % get indices of trials that belong to this ND
            thisAmpInd = find(fictracIInj(i).whichAmp == amps(j));
            % check that there is data
            if ~isempty(thisAmpInd)
                % loop through all FicTrac params
                for k = 1:length(ftParamNames)
                    % FicTrac values
                    fictracIInjMean(i,j).(ftParamNames{k}) = ...
                        mean(fictracIInj(i).(ftParamNames{k})(:,thisAmpInd),2);
                    fictracIInjStdDev(i,j).(ftParamNames{k}) = ...
                        std(fictracIInj(i).(ftParamNames{k})(:,thisAmpInd),[],2);
                    fictracIInjSEM(i,j).(ftParamNames{k}) = ...
                        std(fictracIInj(i).(ftParamNames{k})(:,thisAmpInd),[],2) / ...
                        sqrt(length(thisAmpInd));

                    % change in FicTrac values
                    changeFictracIInjMean(i,j).(ftParamNames{k}) = ...
                        mean(changeFictracIInj(i).(ftParamNames{k})(:,thisAmpInd),2);
                    changeFictracIInjStdDev(i,j).(ftParamNames{k}) = ...
                        std(changeFictracIInj(i).(ftParamNames{k})(:,thisAmpInd),[],2);
                    changeFictracIInjSEM(i,j).(ftParamNames{k}) = ...
                        std(changeFictracIInj(i).(ftParamNames{k})(:,thisAmpInd),[],2) / ...
                        sqrt(length(thisAmpInd));
                end
            end
        end
    end


    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_FicTracOpto.mat'];
    save(saveFileFullName, 'fictracIInj', 'fictracIInjMean', ...
        'fictracIInjStdDev', 'fictracIInjSEM', 'changeFictracIInj', ...
        'changeFictracIInjMean', 'changeFictracIInjStdDev', ...
        'changeFictracIInjSEM', 'trialWindow', ...
        'walkTime', 'minWalkFwd', 'flipLR', 'durs', 'amps', ...
        'trialTimes', '-v7.3');

end