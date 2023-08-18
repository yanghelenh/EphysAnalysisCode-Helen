% extractFicTracOpto_fly.m
%
% Function to extract FicTrac values with optogenetic 
%  stimulation.
% Like extractContStepParamsOpto_fly(), allows filtering valid opto trials
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
%   flipLR - boolean for whether to flip left/right asymmetric vars
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/18/23 - HHY
%
% UPDATED:
%   8/18/23 - HHY
%
function extractFicTracOpto_fly(durs, NDs, trialWindow, walkTime, ...
    minWalkFwd, flipLR, pDataPath, saveFileDir)

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
    for i = 1:length(ftParamNames)
        ftOptoOne.(ftParamNames{i}) = [];
    end
    % to keep track of which ND for each trial
    ftOptoOne.whichND = [];

    % make struct array, 1 for FicTrac values, 1 for change from start
    %  1 struct for each duration
    fictracOpto = repmat(ftOptoOne, numDurs,1);
    changeFictracOpto = repmat(ftOptoOne, numDurs,1);



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

        % check if this pData file has opto, fictracProc structs, if not, 
        %  skip
        if (~any(strcmpi(pDatVarsNames, 'opto')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % load data
        load(pDataFullPath, 'opto', 'fictracProc');

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
                    fictracProc.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(opto.stimCmdDurs(j) == durs);
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
                    % add this trial's FicTrac to running tracker of trials
                    for k = 1:length(ftParamNames)
                        % flip parameter values if needed
                         if (flipLR)
                            % those parameters that need to be value inverted
                            if any(strcmpi(ftParamNames{k}, flipParams))
                                % FicTrac values
                                fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1);

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd) * -1);

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1) - ...
                                    preOptoMean);

                            else
                                % FicTrac values
                                fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preOptoMean);
                            end
                         else
                            fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preOptoMean);
                         end
                    end
                    % add this trial's ND
                    fictracOpto(thisDurInd).whichND = [...
                        fictracOpto(thisDurInd).whichND; ...
                        opto.stimParams.ndFilter];
                    changeFictracOpto(thisDurInd).whichND = [...
                        changeFictracOpto(thisDurInd).whichND; ...
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
                    fictracProc.t, 1, 'last');
                % get end index, based on duration
                thisDurInd = find(opto.stimCmdDurs(j) == durs);
                
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
                                fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1);

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd) * -1);

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd) * -1) - ...
                                    preOptoMean);
                            else
                                % FicTrac values
                                fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd));

                                % change in FicTrac values
                                % get mean before opto
                                preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialPreEndInd));

                                % mean subtracted FicTrac value
                                changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                    changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                    (fictracProc.(ftParamNames{k})(...
                                    thisTrialStartInd:thisTrialEndInd)) - ...
                                    preOptoMean);
                            end
                         else
                            % FicTrac values
                            fictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                fictracOpto(thisDurInd).(ftParamNames{k}), ...
                                fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd));

                            % change in FicTrac values
                            % get mean before opto
                            preOptoMean = mean(fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialPreEndInd));

                            % mean subtracted FicTrac value
                            changeFictracOpto(thisDurInd).(ftParamNames{k}) = cat(2,...
                                changeFictracOpto(thisDurInd).(ftParamNames{k}), ...
                                (fictracProc.(ftParamNames{k})(...
                                thisTrialStartInd:thisTrialEndInd)) - ...
                                preOptoMean);
                         end
                    end
                    % add this trial's ND
                    fictracOpto(thisDurInd).whichND = [...
                        fictracOpto(thisDurInd).whichND; -1];
                    changeFictracOpto(thisDurInd).whichND = [...
                        changeFictracOpto(thisDurInd).whichND; -1];
                end
            end
        end
    end


    % compute means, std dev, SEM
    
    % preallocate
    for i = 1:length(ftParamNames)
        for j = 1:numDurs
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
    
    % loop through all durations
    for i = 1:numDurs
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


    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_FicTracOpto.mat'];
    save(saveFileFullName, 'fictracOpto', 'fictracOptoMean', ...
        'fictracOptoStdDev', 'fictracOptoSEM', 'changeFictracOpto', ...
        'changeFictracOptoMean', 'changeFictracOptoStdDev', ...
        'changeFictracOptoSEM', 'trialWindow', ...
        'walkTime', 'minWalkFwd', 'flipLR', 'durs', 'NDs', ...
        'trialTimes', '-v7.3');

end