% computeLegFictracEphysIinj_1Fly.m
%
% Function that extracts the individual leg traces, FicTrac, and ephys in 
%  response to current injection of different durations and amplitudes for 
%  a single fly.
% User selects all pData for that fly through the file selection GUI
% Saves the data as a struct in path specified in function call
%
% INPUTS:
%   amps - vector of all current injection amplitudes (in pA) to consider
%   durs - vector of all durations of stimulation to consider
%   bwStimDur - scalar value of time between stimulations to consider, in
%       seconds, rep will be stimulation plus this time before and after
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
%
function computeLegFictracEphysIinj_1Fly(amps, durs, bwStimDur, savePath)

    % list of FicTrac variables to process - these are all those in
    %  fictracProc - not the most robust implementation
    fictracVarNames = {'fwdCumPos', 'fwdVel', 'slideCumPos', 'slideVel',...
        'yawAngCumPos','yawAngPosWrap','yawAngVel','yawAngSpd',...
        'totAngSpd','totSpd','xPos','yPos'};
    % list of leg tracking variables to process, from legTrack
    legTrackVarNames = {'srnfLegX', 'srnfLegY', 'legXVel', 'legYVel'};
    % list of ephys variables to process, from ephysData and ephysSpikes
    ephysVarNames = {'scaledVoltage', 'spikeRate', 'medFiltV'};

    NUM_LEG_PTS = 11; % number of tracked points in leg video

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
                any(contains(pDataVarsNames,'ephysSpikes')))

            % load pData
            load(pDataFullPath, 'legTrack','fictracProc', 'iInj', ...
                'ephysData', 'ephysSpikes');

            % loop through all FicTrac variables
            for j = 1:length(fictracVarNames)
                thisVarName = fictracVarNames{j};
                thisVarVal = fictracProc.(thisVarName);

                % replace dropped indices with NaN
                thisVarVal(fictracProc.dropInd) = nan;

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
                    fictracProc.t, iInj, amps, durs, bwStimDur, ...
                    pDataName, norm2StimStart);
            end
            
            % loop through all legTrack variables
            for j = 1:length(legTrackVarNames)
                thisVarName = legTrackVarNames{j};
                thisVarVal = legTrack.(thisVarName);

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
                        thisVarVal(:,k), legTrack.t, iInj, amps, durs, ...
                        bwStimDur, pDataName, false);
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
                    thisVarVal, thisVarT, iInj, amps, durs, bwStimDur, ...
                    pDataName, false);
            end

            % clear this trial's loaded data
            clear fictracProc iInj legTrack ephysData ephysSpikes
        end
    end

    % compute means and std errs for FicTrac vars
    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        [fictracIinj.(thisVarName).means, ...
            fictracIinj.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(fictracIinj.(thisVarName).reps);
    end

    % compute means and std errs for leg tracking vars
    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};

        % loop through all tracked points
        for j = 1:NUM_LEG_PTS
            [legTrackIinj.(thisVarName).means(:,:,j), ...
                legTrackIinj.(thisVarName).stdErrs(:,:,j)] = ...
                computeMeansStdErrsFictracOpto(...
                legTrackIinj.(thisVarName).reps(:,:,j));
        end
    end

    % compute means and std errs for ephys vars
    for i = 1:length(ephysVarNames)
        thisVarName = ephysVarNames{i};
        [ephysIinj.(thisVarName).means, ...
            ephysIinj.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(ephysIinj.(thisVarName).reps);
    end

    % add timing vectors to all data structs
    fictracIinj.durTs = ftDurTs;
    legTrackIinj.durTs = legDurTs;
    ephysIinj.durTs = ephysDurTs;

    % full path to save data
    saveFullPath = [savePath filesep flyName '_avgLegFictracEphysIinj.mat'];

    % save data
    save(saveFullPath, 'fictracIinj','legTrackIinj','ephysIinj', ...
        'amps','durs', 'bwStimDur', '-v7.3');

    fprintf('Saved legFictracEphysIinj for %s!\n', flyName);
end