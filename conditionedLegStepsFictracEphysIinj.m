% conditionedLegStepsFictracEphysIinj.m
%
% Function that takes reps extracted from
%  computeLegStepsFictracEphysIinj_1Fly() and returns only those that 
%  meet the specified criteria: greater than or less than a threshold value
%  on the variable to condition on, during a time window prioer to
%  stimulation
% Returns structs identical to those from 
%  computeLegStepsFictracEphysIinj_1Fly(), with only relevant reps included
%
% INPUTS:
%   datPath - full path to .mat file that's the output of 
%       computeLegStepsFictracEphysIinj_1Fly()
%   condVarParams - struct of parameters for conditioning on
%       name - name of variable to condition on, must be part of ephysIinj,
%           fictracIinj, legStepsIinj, legTrackIinj
%       whichStruct - 'ephysIinj', 'fictracIinj', 'legStepsIinj', or 
%           'legTrackIinj'; specify which struct variable to condition on
%           belongs to
%       legInd - if variable to condition on is part of legStepsIinj or
%           legTrackIinj, then specify which leg to condition on, using 1:6
%           index
%       timeBf - time, in sec, before stimulus start to consider
%           conditioned variable; this is start of time window
%       winDur - duration, in sec, of the window to consider for
%           conditioning
%       whichDir - 'less' or 'greater'; whether to return reps with
%           conditioning variable less than or greater than threshold
%       thresh - threshold, on average value of conditioning variable
%           within the specified time window
%   savePath - full path to where to save output
%
% OUTPUTS:
%   none, but saves file of same name as input file, with _cond appended,
%       containing all the same structs with only the reps that meet the
%       appropriate criteria; also contains condVarParams
%
% CREATED: 8/19/22 - HHY
%
% UPDATED:
%   8/19/22 - HHY
%
function conditionedLegStepsFictracEphysIinj(datPath, condVarParams, ...
    savePath)

    % list of FicTrac variables to process - these are all those in
    %  fictracProc - not the most robust implementation
    fictracVarNames = {'fwdCumPos', 'fwdVel', 'slideCumPos', 'slideVel',...
        'yawAngCumPos','yawAngPosWrap','yawAngVel','yawAngSpd',...
        'totAngSpd','totSpd','xPos','yPos'};
    % list of leg tracking variables to process, from legTrack
    legTrackVarNames = {'srnfLegX', 'srnfLegY', 'legXVel', 'legYVel'};
    % list of ephys variables to process, from ephysData and ephysSpikes
    ephysVarNames = {'scaledVoltage', 'spikeRate', 'medFiltV'};
    % list of leg step variables to process, from legSteps
    legStepsVarNames = {'stepLengths', 'stepXLengths', 'stepYLengths',...
        'stepDirections', 'stepDurations', 'stepSpeeds', 'stepVelX', ...
        'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', 'stepPEPY'};

    NUM_LEG_PTS = 11; % number of tracked points in leg video
    NUM_LEGS = 6; % number of legs
    NUM_PHASES = 2; % swing/stance

    % load in processed data from computeLegStepsFictracEphysIinj_1Fly()
    %  assumes path is correct and all expected data is present and in the
    %  right format
    load(datPath, 'amps', 'bwStimDur', 'durs', 'ephysIinj', ...
        'fictracIinj', 'legStepsIinj', 'legTrackIinj');

    % convert name and whichStruct into actual variable to condition on
    switch condVarParams.whichStruct
        case 'ephysIinj'
            condVarStrct = ephysIinj.(condVarParams.name);
            condVarStrct.durTs = ephysIinj.durTs;
        case 'fictracIinj'
            condVarStrct = fictracIinj.(condVarParams.name);
            condVarStrct.durTs = fictracIinj.durTs;
        case 'legTrackIinj'
            condVarStrctTemp = legTrackIinj.(condVarParams.name);
            
            % get vars for this leg only
            condVarStrct.reps = ...
                condVarStrctTemp.reps(:,:,condVarParams.legInd);
            condVarStrct.repsIinjTimes = ...
                condVarStrctTemp.repsIinjTimes(:,:,condVarParams.legInd);
            condVarStrct.repsPDataNames = ...
                condVarStrctTemp.repsPDataNames(:,:,condVarParams.legInd);
            condVarStrct.durTs = legTrackIinj.durTs;
        case 'legStepsIinj'
            condVarStrctTemp = legStepsIinj.(condVarParams.name);

            % get vars for this leg only
            condVarStrct.reps = ...
                condVarStrctTemp.reps(:,:,condVarParams.legInd);
            condVarStrct.repsIinjTimes = ...
                condVarStrctTemp.repsIinjTimes(:,:,condVarParams.legInd);
            condVarStrct.repsPDataNames = ...
                condVarStrctTemp.repsPDataNames(:,:,condVarParams.legInd);
            condVarStrct.durTs = legStepsIinj.durTs;
        otherwise
            disp('Invalid condVarParams.whichStrct. Ending');
            return;
    end

    % get timing 

    % initialize struct containing iInj times and pData names for every
    % valid rep
    validReps.respIinjTimes = cell(length(amps),length(durs));
    validReps.respPDataNames = cell(length(amps),length(durs));

    % loop through every rep in condVarStrct, evaluate whether it meets
    %  criteria
    for i = 1:length(amps)
        for j = 1:length(durs)
            % matrix of these reps, # reps x rep duration
            theseReps = condVarStrct.reps{i,j};
            % number of reps
            numReps = size(theseReps,1);

            % loop through all reps
            for k = 1:numReps
                % get value of conditioning variable, within time window
                % specified, for this rep
                thisRep = theseReps(k,:); % vector, this rep
                % vector, this rep's timing
                thisRepTime = condVarStrct.durTs{j};

                % get value of cond var within time window
                startInd = find(thisRepTime >= ...
                    (-1 * condVarParams.timeBf),1,'first');
                endInd = find(thisRepTime <= (-1 * condVarParams.timeBf + ...
                    condVarParams.winDur),1,'last');

                avgWinVal = mean(thisRep(startInd:endInd));

                % check if this rep meets criteria (less than or greater
                %  than thresh)
                % switch based on whether greater than or less than thresh
                switch condVarParams.whichDir
                    case 'less'
                        % if meets criteria, append trial info to valid
                        %  strct
                        if (avgWinVal <= condVarParams.thresh)
                            validReps.respIinjTimes{i,j} = [...
                                validReps.respIinjTimes{i,j}; ...
                                condVarStrct.respIinjTimes{i,j}(k)];
                            validReps.respPDataNames{i,j}{end + 1} = ...
                                condVarStrct.respPDataNames{i,j}{k};
                        end
                    case 'greater'
                        % if meets criteria, append trial to valid strct
                        if (avgWinVal >= condVarParams.thresh)
                            validReps.respIinjTimes{i,j} = [...
                                validReps.respIinjTimes{i,j}; ...
                                condVarStrct.repsIinjTimes{i,j}(k)];
                            validReps.respPDataNames{i,j}{end + 1} = ...
                                condVarStrct.repsPDataNames{i,j}{k};
                        end
                end
            end
        end
    end

    % from info on valid reps based on conditioning variable, extract valid
    %  reps for all other variables 

    % preallocate
    fictracFieldStrct.reps = cell(length(amps), length(durs));
    fictracFieldStrct.repsIinjTimes = cell(length(amps), length(durs));
    fictracFieldStrct.repsPDataNames = cell(length(amps), length(durs));
    fictracFieldStrct.means = cell(length(amps), length(durs));
    fictracFieldStrct.stdErrs = cell(length(amps), length(durs));

    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        fictracCond.(thisVarName) = fictracFieldStrct;
    end

    % loop through all fields of FicTrac struct
    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};

        % loop through all stimulation amplitudes
        for j = 1:length(amps)
            % loop through all stimulation durations
            for k = 1:length(durs)
                % matrix of these reps, # reps x rep duration
                theseReps = fictracIinj.(thisVarName).reps{j,k};
                % number of reps
                numReps = size(theseReps,1);

                % loop through all reps
                for l = 1:numReps
                    % get stim start time for this rep
                    thisRepIinjTime = fictracIinj.(thisVarName).repsIinjTimes{j,k}(l);
                    % get pData name for this rep
                    thisRepPDatName = fictracIinj.(thisVarName).repsPDataNames{j,k}{l};

                    % check if this rep is a valid one
                    % check if this stim start time is in valid list
                    iInjTimeInd = find(thisRepIinjTime == ...
                        validReps.respIinjTimes{j,k});
                    if (~isempty(iInjTimeInd))
                        % get pData names for these iInj times
                        thesePDat = validReps.respPDataNames{j,k}(iInjTimeInd);

                        % check if this pDat name is in valid list
                        % if yes, add to reps
                        if (any(strcmp(thisRepPDatName, thesePDat)))
                            fictracCond.(thisVarName).reps{j,k} = [...
                                fictracCond.(thisVarName).reps{j,k}; ...
                                theseReps(l,:)];

                            % also add pData name and iInj time
                            fictracCond.(thisVarName).repsIinjTimes{j,k} = [...
                                fictracCond.(thisVarName).repsIinjTimes{j,k};...
                                thisRepIinjTime];
                            fictracCond.(thisVarName).repsPDataNames{j,k}{end + 1} = ...
                                thisRepPDatName;
                        end
                    end
                end
            end
        end
    end

    % loop through all fields of legTrack struct
    
    % preallocate
    legFieldStrct.reps = cell(length(amps),length(durs),NUM_LEG_PTS);
    legFieldStrct.repsIinjTimes = cell(length(amps),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.repsPDataNames = cell(length(amps),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.means = cell(length(amps),length(durs),NUM_LEG_PTS);
    legFieldStrct.stdErrs = cell(length(amps),length(durs),NUM_LEG_PTS);

    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};
        legTrackCond.(thisVarName) = legFieldStrct;
    end

        % loop through all fields of legTrack struct
    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};

        % loop through all stimulation amplitudes
        for j = 1:length(amps)
            % loop through all stimulation durations
            for k = 1:length(durs)
                % loop through all legs
                for m = 1:NUM_LEG_PTS
                    % matrix of these reps, # reps x rep duration
                    theseReps = legTrackIinj.(thisVarName).reps{j,k,m};
                    % number of reps
                    numReps = size(theseReps,1);
    
                    % loop through all reps
                    for l = 1:numReps
                        % get stim start time for this rep
                        thisRepIinjTime = legTrackIinj.(thisVarName).repsIinjTimes{j,k,m}(l);
                        % get pData name for this rep
                        thisRepPDatName = legTrackIinj.(thisVarName).repsPDataNames{j,k,m}{l};
    
                        % check if this rep is a valid one
                        % check if this stim start time is in valid list
                        iInjTimeInd = find(thisRepIinjTime == ...
                            validReps.respIinjTimes{j,k});
                        if (~isempty(iInjTimeInd))
                            % get pData names for these iInj times
                            thesePDat = validReps.respPDataNames{j,k}(iInjTimeInd);
    
                            % check if this pDat name is in valid list
                            % if yes, add to reps
                            if (any(strcmp(thisRepPDatName, thesePDat)))
                                legTrackCond.(thisVarName).reps{j,k,m} = [...
                                    legTrackCond.(thisVarName).reps{j,k,m}; ...
                                    theseReps(l,:)];
    
                                % also add pData name and iInj time
                                legTrackCond.(thisVarName).repsIinjTimes{j,k,m} = [...
                                    legTrackCond.(thisVarName).repsIinjTimes{j,k,m};...
                                    thisRepIinjTime];
                                legTrackCond.(thisVarName).repsPDataNames{j,k,m}{end + 1} = ...
                                    thisRepPDatName;
                            end
                        end
                    end
                end
            end
        end
    end


    % loop through all fields of ephysVar struct

    % preallocate
    ephysFieldStrct.reps = cell(length(amps), length(durs));
    ephysFieldStrct.repsIinjTimes = cell(length(amps), length(durs));
    ephysFieldStrct.repsPDataNames = cell(length(amps), length(durs));
    ephysFieldStrct.means = cell(length(amps), length(durs));
    ephysFieldStrct.stdErrs = cell(length(amps), length(durs));

    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};
        ephysCond.(thisVarName) = ephysFieldStrct;
    end

    % loop through all fields of ephys struct
    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};

        % loop through all stimulation amplitudes
        for j = 1:length(amps)
            % loop through all stimulation durations
            for k = 1:length(durs)
                % matrix of these reps, # reps x rep duration
                theseReps = ephysIinj.(thisVarName).reps{j,k};
                % number of reps
                numReps = size(theseReps,1);

                % loop through all reps
                for l = 1:numReps
                    % get stim start time for this rep
                    thisRepIinjTime = ephysIinj.(thisVarName).repsIinjTimes{j,k}(l);
                    % get pData name for this rep
                    thisRepPDatName = ephysIinj.(thisVarName).repsPDataNames{j,k}{l};

                    % check if this rep is a valid one
                    % check if this stim start time is in valid list
                    iInjTimeInd = find(thisRepIinjTime == ...
                        validReps.respIinjTimes{j,k});
                    if (~isempty(iInjTimeInd))
                        % get pData names for these iInj times
                        thesePDat = validReps.respPDataNames{j,k}(iInjTimeInd);

                        % check if this pDat name is in valid list
                        % if yes, add to reps
                        if (any(strcmp(thisRepPDatName, thesePDat)))
                            ephysCond.(thisVarName).reps{j,k} = [...
                                ephysCond.(thisVarName).reps{j,k}; ...
                                theseReps(l,:)];

                            % also add pData name and iInj time
                            ephysCond.(thisVarName).repsIinjTimes{j,k} = [...
                                ephysCond.(thisVarName).repsIinjTimes{j,k};...
                                thisRepIinjTime];
                            ephysCond.(thisVarName).repsPDataNames{j,k}{end + 1} = ...
                                thisRepPDatName;
                        end
                    end
                end
            end
        end
    end

    % loop through all fields of legSteps struct

    legStepsFieldStruct.reps = cell(length(amps),length(durs),NUM_LEGS,...
        NUM_PHASES);
    legStepsFieldStruct.repsIinjTimes = cell(length(amps),length(durs),...
        NUM_LEGS,NUM_PHASES);
    legStepsFieldStruct.repsPDataNames = cell(length(amps),length(durs),...
        NUM_LEGS,NUM_PHASES);
    legStepsFieldStrct.means = cell(length(amps),length(durs),...
        NUM_LEGS, NUM_PHASES);
    legStepsFieldStrct.stdErrs = cell(length(amps),length(durs),...
        NUM_LEGS,NUM_PHASES);

    for i = 1:length(legStepsVarNames)
        thisVarName = legStepsVarNames{i};
        legStepsCond.(thisVarName) = legStepsFieldStruct;
    end

        % loop through all fields of legSteps struct
    for i = 1:length(legStepsVarNames)
        thisVarName = legStepsVarNames{i};

        % loop through all stimulation amplitudes
        for j = 1:length(amps)
            % loop through all stimulation durations
            for k = 1:length(durs)
                % loop through all legs
                for m = 1:NUM_LEGS
                    % loop through all phases
                    for n = 1:NUM_PHASES
                        % matrix of these reps, # reps x rep duration
                        theseReps = legStepsIinj.(thisVarName).reps{j,k,m,n};
                        % number of reps
                        numReps = size(theseReps,1);
        
                        % loop through all reps
                        for l = 1:numReps
                            % get stim start time for this rep
                            thisRepIinjTime = legStepsIinj.(thisVarName).repsIinjTimes{j,k,m,n}(l);
                            % get pData name for this rep
                            thisRepPDatName = legStepsIinj.(thisVarName).repsPDataNames{j,k,m,n}{l};
        
                            % check if this rep is a valid one
                            % check if this stim start time is in valid list
                            iInjTimeInd = find(thisRepIinjTime == ...
                                validReps.respIinjTimes{j,k});
                            if (~isempty(iInjTimeInd))
                                % get pData names for these iInj times
                                thesePDat = validReps.respPDataNames{j,k}(iInjTimeInd);
        
                                % check if this pDat name is in valid list
                                % if yes, add to reps
                                if (any(strcmp(thisRepPDatName, thesePDat)))
                                    legStepsCond.(thisVarName).reps{j,k,m,n} = [...
                                        legStepsCond.(thisVarName).reps{j,k,m,n}; ...
                                        theseReps(l,:)];
        
                                    % also add pData name and iInj time
                                    legStepsCond.(thisVarName).repsIinjTimes{j,k,m,n} = [...
                                        legStepsCond.(thisVarName).repsIinjTimes{j,k,m,n};...
                                        thisRepIinjTime];
                                    legStepsCond.(thisVarName).repsPDataNames{j,k,m}{end + 1} = ...
                                        thisRepPDatName;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % compute means and std errs for FicTrac vars
    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        [fictracCond.(thisVarName).means, ...
            fictracCond.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(fictracCond.(thisVarName).reps);
    end

    % compute means and std errs for leg tracking vars
    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};

        % loop through all tracked points
        for j = 1:NUM_LEG_PTS
            [legTrackCond.(thisVarName).means(:,:,j), ...
                legTrackCond.(thisVarName).stdErrs(:,:,j)] = ...
                computeMeansStdErrsFictracOpto(...
                legTrackCond.(thisVarName).reps(:,:,j));
        end
    end

    % compute means and std errs for ephys vars
    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};
        [ephysCond.(thisVarName).means, ...
            ephysCond.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(ephysCond.(thisVarName).reps);
    end

    % compute means and std errs for leg step vars
    for i = 1:length(legStepsVarNames)
        thisVarName = legStepsVarNames{i};

        % loop through all legs
        for j = 1:NUM_LEGS
            % loop through swing/stance
            for k = 1:NUM_PHASES
                [legStepsCond.(thisVarName).means(:,:,j,k),...
                    legStepsCond.(thisVarName).stdErrs(:,:,j,k)] = ...
                    computeMeansStdErrsFictracOpto(...
                    legStepsCond.(thisVarName).reps(:,:,j,k));
            end
        end
    end

    % add durTs to structs
    fictracCond.durTs = fictracIinj.durTs;
    ephysCond.durTs = ephysIinj.durTs;
    legTrackCond.durTs = legTrackIinj.durTs;
    legStepsCond.durTs = legStepsIinj.durTs;

    % get file name
    nameStartInd = find(datPath==filesep, 1,'last') + 1;
    % just file name, without .mat
    datFileName = datPath(nameStartInd:(end-4));
    
    % full path to save data
    saveFullPath = [savePath filesep datFileName '_cond.mat'];

    % save data
    save(saveFullPath, 'fictracCond','legTrackCond','ephysCond', ...
        'legStepsCond', 'amps','durs', 'bwStimDur', ...
        'condVarParams', '-v7.3');

    fprintf('Saved legFictracEphysIinj for %s!\n', datFileName);
end
