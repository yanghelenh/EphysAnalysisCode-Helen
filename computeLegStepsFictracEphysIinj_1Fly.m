% computeLegStepsFictracEphysIinj_1Fly.m
%
% Function that extracts the individual leg traces, FicTrac, and ephys in 
%  response to current injection of different durations and amplitudes for 
%  a single fly.
% User selects all pData for that fly through the file selection GUI
% Saves the data as a struct in path specified in function call
% Has option to shift apparent timing of current injection, to extract
%  control time periods
%
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   bwStimDur - scalar value of time between stimulations to consider, in
%       secoamps, rep will be stimulation plus this time before and after
%   iInjShiftTime - in sec, how much to shift actual current injection start
%       times (neg before, pos after), abs val must be less than bwStimDur;
%       0 is actual; other values as control comparison for +/- current
%       injection
%   savePath - path to folder in which to save data
%
% OUTPUTS:
%   none, but saves struct to file of name date_fly_legFictracEphysIinj.mat 
%       into savePath
%
% CREATED: 6/29/22 - HHY
%
% UPDATED:
%   6/29/22 - HHY
%   9/22/22 - HHY - add iInjShiftTime input for extracting control 
%       stimulations
%   7/19/23 - HHY - fix circular stats
%
function computeLegStepsFictracEphysIinj_1Fly(amps, durs, bwStimDur, ...
    iInjShiftTime, savePath)

    legIDs.ind = 1:6; % indicies into raw position matricies for legs
    legIDs.names = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};

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
    % all step parameters that are circular variables - need to use
    %  circular stats - 6/14/23 - HHY
    circStepParams = {'stepDirections'};

    NUM_LEG_PTS = 11; % number of tracked points in leg video
    NUM_LEGS = 6; % number of legs
    NUM_PHASES = 2; % swing/stance

    STEP_INTERP_FR = 250; % step interpolation frame rate

    % prompt user to select pData files
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate struct for FicTrac data
    % for each variable, struct for each with reps, means, and std errors
    % reps as #amps x #durs cell array, with each element being a matrix of
    %  reps x rep length 
    % repsIinjTimes as #amps x #durs cell array, with each element being a
    %  vector of # reps length
    % repsPDataNames as #amps x #durs cell array, with each element being a
    %  cell array of # reps length
    % means and stdErrs as same size cell array, each element as vector 
    fictracFieldStrct.reps = cell(length(amps), length(durs));
    fictracFieldStrct.repsIinjTimes = cell(length(amps), length(durs));
    fictracFieldStrct.repsPDataNames = cell(length(amps), length(durs));
    fictracFieldStrct.means = cell(length(amps), length(durs));
    fictracFieldStrct.stdErrs = cell(length(amps), length(durs));

    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        fictracIinj.(thisVarName) = fictracFieldStrct;
    end

    % preallocate struct for leg data
    % for each variable, struct for each rep, mean, and std error
    % reps as #amps x #durs x #tracked points cell array, with each element
    %  being a matrix of reps x rep length
    % repsInjTimes as #amps x #durs x #tracked points cell array, with each
    %  element being a vector of # reps length, tracking opto stim start
    %  time corresponding to each rep
    % repsPDataNames as #amps x #durs x #tracked points cell array, with
    %  each element being a cell array of # reps length, tracking the name
    %  of the pData file the rep came from
    % means and stdErrs as same size cell array, each element as vector
    legFieldStrct.reps = cell(length(amps),length(durs),NUM_LEG_PTS);
    legFieldStrct.repsIinjTimes = cell(length(amps),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.repsPDataNames = cell(length(amps),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.means = cell(length(amps),length(durs),NUM_LEG_PTS);
    legFieldStrct.stdErrs = cell(length(amps),length(durs),NUM_LEG_PTS);

    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};
        legTrackIinj.(thisVarName) = legFieldStrct;
    end

    % preallocate struct for ephys data
    % for each variable, struct for each with reps, means, and std errors
    % reps as #amps x #durs cell array, with each element being a matrix of
    %  reps x rep length 
    % repsOptoTimes as #amps x #durs cell array, with each element being a
    %  vector of # reps length
    % repsPDataNames as #amps x #durs cell array, with each element being a
    %  cell array of # reps length
    % means and stdErrs as same size cell array, each element as vector 
    ephysFieldStrct.reps = cell(length(amps), length(durs));
    ephysFieldStrct.repsIinjTimes = cell(length(amps), length(durs));
    ephysFieldStrct.repsPDataNames = cell(length(amps), length(durs));
    ephysFieldStrct.means = cell(length(amps), length(durs));
    ephysFieldStrct.stdErrs = cell(length(amps), length(durs));

    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};
        ephysIinj.(thisVarName) = ephysFieldStrct;
    end

    % preallocate struct for leg step data
    % for each variable, struct for each rep, mean, and std error
    % reps as #amps x #durs x 6 (# legs) x 2 (swing/stance) cell array,
    %  with each element being a matrix of reps x rep length
    % repsOptoTimes as #amps x #durs x 6 x 2 cell array, with each
    %  element being a vector of # reps length, tracking opto stim start
    %  time corresponding to each rep
    % repsPDataNames as #amps x #durs x 6 x 2 cell array, with
    %  each element being a cell array of # reps length, tracking the name
    %  of the pData file the rep came from
    % means and stdErrs as same size cell array, each element as vector
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
        legStepsIinj.(thisVarName) = legStepsFieldStruct;
    end

    
    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        % save fly name as first pDataName's date and fly (12 characters)
        if (i == 1)
            flyName = pDataName(1:12);
        end
        
        pDataFullPath = [pDataPath pDataName];
        
        % get variables in pData file
        pDataMatObj = matfile(pDataFullPath);
        pDataVarsStrct = whos(pDataMatObj);
        pDataVars = struct2cell(pDataVarsStrct);
        pDataVarsNames = pDataVars(1,:); % cell array of names
        
        % check that pData has legTrack, FicTrac, ephys, and Iinj data, 
        %  otherwise, skip
        if (any(contains(pDataVarsNames, 'legTrack')) && ...
                any(contains(pDataVarsNames, 'iInj')) && ...
                any(contains(pDataVarsNames,'fictracProc')) && ...
                any(contains(pDataVarsNames,'ephysData')) && ...
                any(contains(pDataVarsNames,'ephysSpikes')) && ...
                any(contains(pDataVarsNames,'moveNotMove')) && ...
                any(contains(pDataVarsNames,'legSteps')))

            % load pData
            load(pDataFullPath, 'legTrack','fictracProc', 'iInj', ...
                'ephysData', 'ephysSpikes','moveNotMove','legSteps');

            % get resampling time for legStep
            legStepSampTime = legTrack.t(1):(1/STEP_INTERP_FR):legTrack.t(end);

            % convert iInj to shifted iInj - time offset for stim times,
            %  for ctrl comparisons
            iInjShifted = shiftIinjStartTime(iInj, iInjShiftTime);

            % loop through all FicTrac variables
            for j = 1:length(fictracVarNames)
                thisVarName = fictracVarNames{j};
                thisVarVal = fictracProc.(thisVarName);

                % replace dropped indices with NaN
                thisVarVal(fictracProc.dropInd) = nan;

                % replace not moving indices with NaN
                thisVarVal(moveNotMove.ftNotMoveInd) = nan;

                % if this variable is a position variable, normalize to
                %  value at stimulation start
                if (contains(thisVarName,'pos','IgnoreCase',true))
                    norm2StimStart = true;
                else
                    norm2StimStart = false;
                end
                
                % get all reps for this pData for this FicTrac variable
                [fictracIinj.(thisVarName).reps, ...
                    fictracIinj.(thisVarName).repsIinjTimes, ...
                    fictracIinj.(thisVarName).repsPDataNames,...
                    ftDurTs] = ...
                    extractTrialsIinj(...
                    fictracIinj.(thisVarName).reps, ...
                    fictracIinj.(thisVarName).repsIinjTimes, ...
                    fictracIinj.(thisVarName).repsPDataNames,...
                    thisVarVal, ...
                    fictracProc.t, iInjShifted, amps, durs, bwStimDur, ...
                    pDataName, norm2StimStart);
            end
            
            % loop through all legTrack variables
            for j = 1:length(legTrackVarNames)
                thisVarName = legTrackVarNames{j};
                thisVarVal = legTrack.(thisVarName);

                % replace not moving indices with NaN
                thisVarVal(moveNotMove.legNotMoveInd,:) = nan;

                % loop through all tracked points, get reps
                for k = 1:NUM_LEG_PTS
                    [legTrackIinj.(thisVarName).reps(:,:,k), ...
                        legTrackIinj.(thisVarName).repsIinjTimes(:,:,k),...
                        legTrackIinj.(thisVarName).repsPDataNames(:,:,k),...
                        legDurTs] = ...
                        extractTrialsIinj(...
                        legTrackIinj.(thisVarName).reps(:,:,k), ...
                        legTrackIinj.(thisVarName).repsIinjTimes(:,:,k),...
                        legTrackIinj.(thisVarName).repsPDataNames(:,:,k),...
                        thisVarVal(:,k), legTrack.t, iInjShifted, amps, ...
                        durs, bwStimDur, pDataName, false);
                end
            end

            % loop through all ephys variables
            for l = 1:length(ephysVarNames)
                thisVarName = ephysVarNames{l};

                % since ephys data is distributed between ephysData and
                %  ephysSpikes, make sure to grab the variable values
                %  correctly
                switch thisVarName
                    case 'scaledVoltage'
                        thisVarVal = ephysData.(thisVarName);
                        thisVarT = ephysData.t;
                    otherwise
                        thisVarVal = ephysSpikes.(thisVarName);
                        thisVarT = ephysSpikes.t;
                end

                % get all reps for this pData for this ephys variable
                [ephysIinj.(thisVarName).reps, ...
                    ephysIinj.(thisVarName).repsIinjTimes, ...
                    ephysIinj.(thisVarName).repsPDataNames,...
                    ephysDurTs] = ...
                    extractTrialsIinj(...
                    ephysIinj.(thisVarName).reps, ...
                    ephysIinj.(thisVarName).repsIinjTimes, ...
                    ephysIinj.(thisVarName).repsPDataNames,...
                    thisVarVal, thisVarT, iInjShifted, amps, durs, ...
                    bwStimDur, pDataName, false);
            end

            % loop through all legSteps variables
            for m = 1:length(legStepsVarNames)
                thisVarName = legStepsVarNames{m};
                
                % loop through all legs
                for n = 1:NUM_LEGS
                    whichLeg = legIDs.names{n};
                    whichLegInd = legIDs.ind(n);

                    % loop through both phases (swing/stance)
                    for o = 1:NUM_PHASES

                        switch o
                            case 1
                                whichPhase = 'stance';
                            case 2
                                whichPhase = 'swing';
                        end
                        
                        % get step values over time, constant during step
                        %  duration
                        thisVarVal = stepParam2Vector(legSteps, ...
                            legStepSampTime, whichLeg, thisVarName, ...
                            whichPhase);

                        % get reps
                        [legStepsIinj.(thisVarName).reps(:,:,whichLegInd,o),...
                            legStepsIinj.(thisVarName).repsIinjTimes(:,:,whichLegInd,o),...
                            legStepsIinj.(thisVarName).repsPDataNames(:,:,whichLegInd,o),...
                            legStepsDurTs] = extractTrialsIinj(...
                            legStepsIinj.(thisVarName).reps(:,:,whichLegInd,o), ...
                            legStepsIinj.(thisVarName).repsIinjTimes(:,:,whichLegInd,o),...
                            legStepsIinj.(thisVarName).repsPDataNames(:,:,whichLegInd,o),...
                            thisVarVal, legStepSampTime, iInjShifted, ...
                            amps, durs, bwStimDur, pDataName, false);
                    end
                end
            end

            % clear this trial's loaded data
            clear fictracProc iInj legTrack ephysData ephysSpikes legSteps moveNotMove
        end
    end

    % compute means and std errs for FicTrac vars
    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        [fictracIinj.(thisVarName).means, ...
            fictracIinj.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(fictracIinj.(thisVarName).reps,...
            false);
    end

    % compute means and std errs for leg tracking vars
    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};

        % loop through all tracked points
        for j = 1:NUM_LEG_PTS
            [legTrackIinj.(thisVarName).means(:,:,j), ...
                legTrackIinj.(thisVarName).stdErrs(:,:,j)] = ...
                computeMeansStdErrsFictracOpto(...
                legTrackIinj.(thisVarName).reps(:,:,j), false);
        end
    end

    % compute means and std errs for ephys vars
    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};
        [ephysIinj.(thisVarName).means, ...
            ephysIinj.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(ephysIinj.(thisVarName).reps, ...
            false);
    end

    % compute means and std errs for leg step vars
    for i = 1:length(legStepsVarNames)
        thisVarName = legStepsVarNames{i};

        % loop through all legs
        for j = 1:NUM_LEGS
            % loop through swing/stance
            for k = 1:NUM_PHASES
                if (strcmpi(thisVarName, circStepParams)) % circular
                    [legStepsIinj.(thisVarName).means(:,:,j,k),...
                        legStepsIinj.(thisVarName).stdErrs(:,:,j,k)] = ...
                        computeMeansStdErrsFictracOpto(...
                        legStepsIinj.(thisVarName).reps(:,:,j,k), true);
                else
                    [legStepsIinj.(thisVarName).means(:,:,j,k),...
                        legStepsIinj.(thisVarName).stdErrs(:,:,j,k)] = ...
                        computeMeansStdErrsFictracOpto(...
                        legStepsIinj.(thisVarName).reps(:,:,j,k), false);
                end
            end
        end
    end

    % add timing vectors to all data structs
    fictracIinj.durTs = ftDurTs;
    legTrackIinj.durTs = legDurTs;
    ephysIinj.durTs = ephysDurTs;
    legStepsIinj.durTs = legStepsDurTs;

    % full path to save data, include shift time, after s
    saveFullPath = sprintf('%s%s%s_avgLegStepsFictracEphysIinj_s%d.mat',...
        savePath, filesep, flyName, iInjShiftTime);

    % save data
    save(saveFullPath, 'fictracIinj','legTrackIinj','ephysIinj', ...
        'legStepsIinj', 'amps','durs', 'bwStimDur', '-v7.3');

    fprintf('Saved legFictracEphysIinj for %s!\n', flyName);
end