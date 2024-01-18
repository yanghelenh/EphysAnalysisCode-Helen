% extractFicTracOptoCond_fly.m
%
% Function to extract FicTrac values with optogenetic 
%  stimulation.
% Adaptation of extractFicTracOpto_fly with more flexible conditioning on
%  FicTrac behavior (specify parameter, condition, time relative to opto
%  stim). Note that conditioning on minimum forward walking speed
%  through this instead
% Select all pData files for 1 fly through GUI
% Saves output file with name defined by first pData file (without trial #)
% 
% INPUTS:
%   durs - vector of all durations of stimulation to consider
%   NDs - vector of all NDs to consider
%   trialWindow - length 2 vector where 1st element is time before opto
%       starts to consider as trial and 2nd element is time after opto stim
%       ends to consider as part of the trial
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
%   flipLR - boolean for whether to flip left/right asymmetric vars
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%   
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 12/14/23 - HHY
%
% UPDATED:
%   12/15/23 - HHY
%
function extractFicTracOptoCond_fly(durs, NDs, trialWindow, cond,...
    flipLR, pDataPath, saveFileDir)

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


        % loop through each opto stim trial, check if it meets all cond
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

            % if this trial is valid (meets all cond criteria), record 
            %  trial info
            if (meetsCond)
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
        % check if conditions are met
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


            % if this trial is valid (meets conditions), record trial
            %  info
            if (meetsCond)
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
        'cond', 'flipLR', 'durs', 'NDs', ...
        'trialTimes', '-v7.3');

end