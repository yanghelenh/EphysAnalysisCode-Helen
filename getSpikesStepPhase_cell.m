% getSpikesStepPhase_cell.m
%
% Function to get timing of spikes relative to leg step phase
% Bins spikes by phase during step cycle, then normalizes by time spent
%  during each phase. Specify number of bins of phase, delay between
%  behavior and ephys. Bins separately for each leg
% Can condition on behavioral parameters
% Operates on multiple pData files from 1 cell. Specify through input or
%  select through GUI
%
% INPUTS:
%   cond - struct of conditions that time points have to meet to be
%     included. Multiple conditions are AND. [] if no conditions
%       whichParam - cell array (even if 1 element) on which fictracProc or
%           legStepsCont fields to condition on. Use 'stepFwdBool' for
%           boolean of whether step is moving backwards or forwards. 
%       legs - cell array for which leg for legStepsCont and 'stepFwdBool'.
%         [] for FicTrac elements
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam  
%   tDelay - time offset b/w behavior and ephys. Negative for ephys before
%     behavior. Can be more than 1 delay
%   numPhaseBins - number of bins to divide leg phase into
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataPath - full path to pData files
%   pDataFNames - string or cell array of strings for pData files to operate
%     on. If [], then select pData files through GUI. One fly only 
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves file
%
% CREATED: 9/2/23 - HHY
%
% UPDATED:
%   9/2/23 - HHY
%
function getSpikesStepPhase_cell(cond, tDelay, numPhaseBins, ...
    postStimExclDur, pDataPath, pDataFNames, saveFileName, saveFileDir)

    % preallocate matrices for binning spikes and phase
    numDelays = length(tDelay);
    % 6 for number of legs
    spikeCounts = zeros(numDelays, numPhaseBins, 6);
    phaseCounts = zeros(numDelays, numPhaseBins, 6);

    % get phase bin edges (edges, 1 more than number of bins)
    phaseBinEdges = linspace(0, 360, numPhaseBins + 1);
    phaseBinStarts = phaseBinEdges(1:(end-1));
    phaseBinEnds = phaseBinEdges(2:end);

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

        % check if this pData file has necessary variables, otherwise, skip
        if (~any(strcmpi(pDatVarsNames, 'legPhase')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')) || ...
                ~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'legTrack')))
            continue;
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'legPhase', 'fictracProc', ...
                'legTrack', 'legStepsCont', 'ephysSpikes', 'legSteps', ...
                'iInj');
        else
            load(pDataFullPath, 'legPhase', 'fictracProc', ...
                'legTrack', 'legStepsCont', 'ephysSpikes', 'legSteps');
        end

        % timing vectors
        phaseT = legTrack.t';
        ephysT = ephysSpikes.t;

        % shift phaseT to define bin starts and ends to find where spike
        % falls
        % t is bin start now
        phaseTStarts = phaseT(1:(end-1));
        phaseTEnds = phaseT(2:end);

        % set current injection times to NaN in ephysT
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            % current injection times; only add to end
            iInjStartTimes = iInj.startTimes;
            iInjEndTimes = iInj.endTimes + postStimExclDur;

            % loop through all iInj bouts, set appropriate inds to NaN
            for j = 1:length(iInjStartTimes)
                thisIInjLog = (ephysT > iInjStartTimes(j)) & ...
                    (ephysT < iInjEndTimes(j));

                ephysT(thisIInjLog) = nan;       
            end
        end

        % get spike times
        spikeTimes = ephysT(ephysSpikes.startInd);
        % remove NaNs (which removes spikes during current inj)
        spikeTimes(isnan(spikeTimes)) = [];

        % for leg phase, set points that don't meet behavioral conditions
        %  to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(phaseT));
            for j = 1:length(cond.whichParam)
                % check which data structures conditions belong to
                % check first for stepFwdBool, as it's not a field of any
                %  struct
                if (strcmpi(cond.whichParam{j}, 'stepFwdBool'))
                    thisFwdLog = getStepFwdLogical(legSteps, phaseT,...
                        cond.legs{j});
                    if ~(eval(cond.cond{j})) % if target is false, invert
                        thisLog = ~thisFwdLog;
                    else
                        thisLog = thisFwdLog;
                    end
                % if FicTrac    
                elseif isfield(fictracProc,cond.whichParam{j})
                    thisFT = interp1(fictracProc.t, ...
                        fictracProc.(cond.whichParam{j}), phaseT, ...
                        'spline');
                    thisLog = eval(['thisFT' cond.cond{j}]); 
                % if step parameter    
                elseif isfield(legStepsCont,cond.whichParam{j})
                    thisLegInd = legSteps.legIDs.ind(strcmpi(legSteps.legIDs.names, ...
                        cond.legs{j}));
                    thisStep = interp1(legStepsCont.t, ...
                        legStepsCont.(cond.whichParam{j})(:,thisLegInd), ...
                        phaseT, 'spline');
                    thisLog = eval(['thisStep' cond.cond{j}]);
                end
                % update logical for all conditions
                condLog = condLog & thisLog;
            end
            % invert condition logical
            condFalseLog = ~condLog;
            % set every time when condition logcial false to nan
            legPhase(condFalseLog,:) = nan;
        end


        % loop through all delays
        for j = 1:numDelays
            % adjust spike timing
            thisSpikeTimes = spikeTimes - tDelay(j);

            % loop through all legs
            for k = 1:size(legPhase,2)

                % loop through all spikes, add to counter
                for l = 1:length(thisSpikeTimes)
                    % get index into legPhase 
                    phaseInd = find((thisSpikeTimes(l) >= phaseTStarts) & ...
                        (thisSpikeTimes(l) < phaseTEnds));
                    if ~isempty(phaseInd)
                        % the phase where this spike falls
                        thisSpikePhase = legPhase(phaseInd,k);

                        % if this spike falls during a valid time
                        if ~isnan(thisSpikePhase)
                            % deal with edge case of phase = 360
                            if (thisSpikePhase == 360)
                                thisBinInd = numPhaseBins;
                            else
                            % get phase bin this spike belongs to
                                thisBinInd = find((thisSpikePhase >= phaseBinStarts) & ...
                                    (thisSpikePhase < phaseBinEnds));
                            end

                            if ~isempty(thisBinInd)
                            % update counter
                                spikeCounts(j,thisBinInd, k) = ...
                                    spikeCounts(j,thisBinInd, k) + 1;
                            end
                        end
                    end
                end

                % get this leg phase without nans
                thisLegPhase = legPhase(:,k);
                noNanPhase = thisLegPhase(~isnan(thisLegPhase));

                % loop through all phases
                for l = 1:length(noNanPhase)
                    thisPhase = noNanPhase(l);

                    % deal with edge case of phase = 360
                    if (thisPhase == 360)
                        thisBinInd = numPhaseBins;
                    else
                        thisBinInd = find((thisPhase >= phaseBinStarts) & ...
                            (thisPhase < phaseBinEnds));
                    end

                    if ~isempty(thisBinInd)
                        % update counter
                        phaseCounts(j,thisBinInd,k) = ...
                            phaseCounts(j,thisBinInd,k) + 1;
                    end
                end
            end
        end
    end

    % turn counts into normalized num spikes per phase bin (spike rate)
    phaseTimes = phaseCounts * median(diff(legTrack.t));

    normPhaseSpikes = spikeCounts ./ phaseTimes;

    % save output file
    saveFileFullName = [saveFileDir filesep saveFileName '.mat'];
    save(saveFileFullName, 'normPhaseSpikes', 'spikeCounts', ...
        'phaseCounts', 'phaseTimes', 'tDelay', 'phaseBinEdges', ...
        'cond', 'postStimExclDur', '-v7.3');

end