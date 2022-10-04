% sortLegStepsByOpto.m
%
% Function that categorizes all leg steps by whether they're during an
%  opto stim bout (and opto stim duration) or not.
% Labels every half step with category index number, and generates key for
%  matching b/w category indices and conditions (different durations).
% Not during opto stim has input specifying how long after opto stim has 
%  ended before considering steps as not during opto stim
% Half steps that start before opto but end during are ignored. Half steps
%  that start during opto but end after opto ends are included. Half steps
%  that start before not opto time or end during opto are ignored.
% Different categories for different stimulation durations. Each pData file
%  only has one ND (stim amplitude)
% Opterates on single pData file, selected by GUI or as input full path.
% 
% INPUTS:
%   durs - vector of all durations of stimulation to consider
%   notOptoTime - time in sec after opto stim turns off to not include 
%       in not opto category
%   pDataFullPath - full path to pData file, optional input ([] to use GUI)
%
% OUTPUTS:
%   none, but saves vector of category labels and key back into same pData
%       file
%
% CREATED: 9/30/22 - HHY
%
% UPDATED:
%   9/30/22 - HHY
%   10/3/22 - HHY - change no stim val to -1 from 0 in key
%
function sortLegStepsByOpto(durs, notOptoTime, pDataFullPath)

    % use GUI to select pData file
    if isempty(pDataFullPath)
        % prompt user to select pData files
        [pDataName, pDataPath] = uigetfile('*.mat', 'Select pData file', ...
            pDataDir(), 'MultiSelect', 'off');
        
        pDataFullPath = [pDataPath pDataName];
    else
        % get pDataName
        fileSepInd = find(pDataFullPath == filesep, 1, 'last');
        % pDataName is everything after final filesep
        pDataName = pDataFullPath((fileSepInd+1):end);
    end
        
    % get variables in pData file
    pDataMatObj = matfile(pDataFullPath);
    pDataVarsStrct = whos(pDataMatObj);
    pDataVars = struct2cell(pDataVarsStrct);
    pDataVarsNames = pDataVars(1,:); % cell array of names
    
    % check that pData has legSteps and opto structs, otherwise, skip
    if (any(contains(pDataVarsNames, 'opto')) && ...
        any(contains(pDataVarsNames,'legSteps')))

        % load pData
        load(pDataFullPath, 'opto', 'legSteps');

        % number of steps (across all 6 legs)
        numSteps = length(legSteps.stepWhichLeg);

        % preallocate matrix for keeping track of opto category for all
        %  steps (num steps, treat each half step
        %  separately)
        stepOptoCat = zeros(size(legSteps.stepLengths));
        % values will be based on key below, NaN for not fitting into any
        %  category (boundary steps)

        % build key for opto categories
        % number of durs plus extra 1 for no opto (-1 value)
        numCat = length(durs) + 1;
        optoCat = ones(1,numCat) * -1;
        % counter index into vector, skip 1 for no opto
        counter = 2; 

        % assign durs to indices
        optoCat(2:end) = durs;

        % loop through all steps
        for i = 1:numSteps
            % loop through 2 half steps
            for j = 1:size(legSteps.stepLengths, 2)
                
                % this half step, start and end times
                startTime = legSteps.stepT(i,j);
                endTime = legSteps.stepT(i,j+1);

                % check if this step falls during opto stim, not opto stim,
                %  or neither
                % first opto stim start time that's later than step
                %  start time
                whichOptoInd = find(opto.stimStartTimes < startTime, ...
                    1, 'last');
                % steps before opto stim starts will return empty
                if ~isempty(whichOptoInd)
                    % end time for this opto step
                    thisOptoEndTime = opto.stimEndTimes(whichOptoInd);

                    % get next opto stim start time, if present;
                    % otherwise, set to infinity, for later comparison
                    if (whichOptoInd < length(opto.stimStartTimes))
                        nextOptoStartTime = opto.stimStartTimes(...
                            whichOptoInd + 1);
                    else
                        nextOptoStartTime = inf;
                    end

                    % check that step start time is within opto step:
                    %  less than opto stim end time
                    % step has to start during opto, but can end after opto
                    %  ends
                    if (startTime < thisOptoEndTime)
                        % assign this step to this type of iInj step (by
                        % amp and dur)
                        % get this opto stim's dur (commanded, not actual)
                        thisDur = opto.stimCmdDurs(whichOptoInd);

                        % convert this dur into index
                        thisIndex = find(thisDur == optoCat);

                        % assign this step to this index
                        stepOptoCat(i,j) = thisIndex;

                    % check if this step starts after opto stim ends, 
                    %  including buffer specified by notOptoTime input
                    % and check that the step ends before the next opto
                    %  stim starts
                    elseif (startTime >= (thisOptoEndTime + notOptoTime)...
                            && (endTime < nextOptoStartTime))
                        % if yes, then this is assigned index 1 (by
                        % default)
                        stepOptoCat(i,j) = 1;

                    % if step is not during opto stim and isn't in the not 
                    %  injection window, assign NaN as index 
                    else
                        stepOptoCat(i,j) = nan;
                    end
                % if this step falls before all opto stim (returns empty), 
                %  could fall in not opto stim category if it ends 
                % before opto stim starts
                else
                    % step ends before opto stim starts, is not opto stim
                    if (endTime < opto.stimStartTimes(1))
                        stepOptoCat(i,j) = 1;
                    % otherwise, not in any category, NaN
                    else
                        stepOptoCat(i,j) = nan;
                    end
                end
            end
        end

        % clear loaded data
        clear opto legSteps

        % generate struct for output
        legStepsByOpto.stepOptoCat = stepOptoCat;
        legStepsByOpto.optoCat = optoCat;
        legStepsByOpto.notOptoTime = notOptoTime;

        % save back into same pData file
        save(pDataFullPath, 'legStepsByOpto' ,'-append');

        fprintf('Saved legStepsByOpto into %s\n', pDataName);
    else
        fprintf('%s does not contain opto and legSteps structs\n', ...
            pDataName);
    end
end