% sortLegStepsByIInj.m
%
% Function that categorizes all leg steps by whether they're during a
%  current injection bout (and by injection amplitude) or not.
% Labels every half step with category index number, and generates key for
%  matching b/w category indices and conditions.
% Not during current injection has input specifying how long after current
%  injection has ended before considering steps as not during current
%  injection.
% Half steps that start before iInj but end during are ignored. Half steps
%  that start during iInj but end after iInj ends are included. Half steps
%  that start before not iInj time or end during iInj are ignored.
% Different categories for different injection amplitudes and durations
% Opterates on single pData file, selected by GUI or specified through
%  input as full path.
% 
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   notIinjTime - time in sec after iInj turns off to not include in not
%       iInj category
%   pDataFullPath - full path to pData file, optional input ([] to use GUI)
%
% OUTPUTS:
%   none, but saves vector of category labels and key back into same pData
%       file
%
% CREATED: 9/28/22 - HHY
%
% UPDATED:
%   9/29/22 - HHY
%
function sortLegStepsByIInj(amps, durs, notIinjTime, pDataFullPath)

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
    
    % check that pData has legSteps and iInj structs, otherwise, skip
    if (any(contains(pDataVarsNames, 'iInj')) && ...
        any(contains(pDataVarsNames,'legSteps')))

        % load pData
        load(pDataFullPath, 'iInj', 'legSteps');

        % number of steps (across all 6 legs)
        numSteps = length(legSteps.stepWhichLeg);

        % preallocate matrix for keeping track of current injection
        %  category for all steps (num steps, treat each half step
        %  separately)
        stepIinjCat = zeros(size(legSteps.stepLengths));
        % values will be based on key below, NaN for not fitting into any
        %  category (boundary steps)

        % build key for I inj categories, 1 vector for amps, 1 vector for
        %  durs
        % extra 1 for no iInj 
        numCat = length(amps) * length(durs) + 1;
        iInjCatAmps = zeros(1,numCat);
        iInjCatDurs = zeros(1,numCat);
        % counter index into vectors, skip 1 for 0,0, no iInj
        counter = 2; 

        % assign amps to indices
        for i = 1:length(amps)
            for j = 1:length(durs)
                iInjCatAmps(counter) = amps(i);
                iInjCatDurs(counter) = durs(j);

                counter = counter + 1;
            end
        end

        % loop through all steps
        for i = 1:numSteps
            % loop through 2 half steps
            for j = 1:size(legSteps.stepLengths, 2)
                
                % this half step, start and end times
                startTime = legSteps.stepT(i,j);
                endTime = legSteps.stepT(i,j+1);

                % check if this step falls during current injection, not
                %  current injection, or neither
                % first current injection start time that's later than step
                %  start time
                whichIinjInd = find(iInj.startTimes < startTime,1,'last');
                % steps before iInj starts will return empty
                if ~isempty(whichIinjInd)
                    % end time for this current injection step
                    thisIinjEndTime = iInj.endTimes(whichIinjInd);

                    % get next current injection start time, if present;
                    % otherwise, set to infinity, for later comparison
                    if (whichIinjInd < length(iInj.startTimes))
                        nextIinjStartTime = iInj.startTimes(...
                            whichIinjInd + 1);
                    else
                        nextIinjStartTime = inf;
                    end

                    % check that step start time is within iInj step:
                    %  less than iInj end time
                    % step has to start during iInj, but can end after iInj
                    %  ends
                    if (startTime < thisIinjEndTime)
                        % assign this step to this type of iInj step (by
                        % amp and dur)
                        % get this iInj's amp and dur
                        thisAmp = iInj.amps(whichIinjInd);
                        thisDur = iInj.durs(whichIinjInd);

                        % convert this amp and dur into index
                        thisIndex = intersect(...
                            find(thisAmp == iInjCatAmps), ...
                            find(thisDur == iInjCatDurs));

                        % assign this step to this index
                        stepIinjCat(i,j) = thisIndex;

                    % check if this step starts after current injection
                    % ends, including buffer specified by notIinjTime input
                    % and check that the step ends before the next current
                    %  injection starts
                    elseif (startTime >= (thisIinjEndTime + notIinjTime)...
                            && (endTime < nextIinjStartTime))
                        % if yes, then this is assigned index 1 (by
                        % default)
                        stepIinjCat(i,j) = 1;

                    % if step is not during iInj step and isn't in the not 
                    %  injection window, assign NaN as index 
                    else
                        stepIinjCat(i,j) = nan;
                    end
                % if this step falls before all iInj (returns empty), 
                %  could fall in not current injection category if it ends 
                % before current injection starts
                else
                    % step ends before current injection starts, is not
                    %  iInj
                    if (endTime < iInj.startTimes(1))
                        stepIinjCat(i,j) = 1;
                    % otherwise, not in any category, NaN
                    else
                        stepIinjCat(i,j) = nan;
                    end
                end
            end
        end

        % clear loaded data
        clear iInj legSteps

        % generate struct for output
        legStepsByIinj.stepIinjCat = stepIinjCat;
        legStepsByIinj.iInjCatAmps = iInjCatAmps;
        legStepsByIinj.iInjCatDurs = iInjCatDurs;
        legStepsByIinj.notIinjTime = notIinjTime;

        % save back into same pData file
        save(pDataFullPath, 'legStepsByIinj' ,'-append');

        fprintf('Saved legStepsByIinj into %s\n', pDataName);
    else
        fprintf('%s does not contain iInj and legSteps structs\n', ...
            pDataName);
    end
end