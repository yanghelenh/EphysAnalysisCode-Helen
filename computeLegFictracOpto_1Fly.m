% computeLegFictracOpto_1Fly.m
%
% Function that extracts the individual leg traces and FicTrac, in response
%  to optogenetic stimulation of different durations and intensities (ND
%  filters) for a single fly.
% User selects all pData for that fly through the file selection GUI
% Saves the data as a struct in path specified in function call
%
% INPUTS:
%   NDs - vector of all NDs to consider
%   durs - vector of all durations of stimulation to consider
%   bwStimDur - scalar value of time between stimulations to consider, in
%       seconds, rep will be stimulation plus this time before and after
%   savePath - path to folder in which to save data
%
% OUTPUTS:
%   none, but saves struct to file of name date_fly_legOpto.mat into
%       savePath
%
% CREATED: 5/12/22 - HHY
%
% UPDATED:
%   5/12/22 - HHY
%   5/19/22 - HHY - add tracking of opto stim start time for each rep
%   5/20/22 - HHY - add tracking of pData name for each rep
%
function computeLegFictracOpto_1Fly(NDs, durs, bwStimDur, savePath)

    % list of FicTrac variables to process - these are all those in
    %  fictracProc - not the most robust implementation
    fictracVarNames = {'fwdCumPos', 'fwdVel', 'slideCumPos', 'slideVel',...
        'yawAngCumPos','yawAngPosWrap','yawAngVel','yawAngSpd',...
        'totAngSpd','totSpd','xPos','yPos'};
    % list of leg tracking variables to process, from legTrack
    legTrackVarNames = {'srnfLegX', 'srnfLegY', 'legXVel', 'legYVel'};

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
    % reps as #NDs x #durs cell array, with each element being a matrix of
    %  reps x rep length 
    % repsOptoTimes as #NDs x #durs cell array, with each element being a
    %  vector of # reps length
    % repsPDataNames as #NDs x #durs cell array, with each element being a
    %  cell array of # reps length
    % means and stdErrs as same size cell array, each element as vector 
    fictracFieldStrct.reps = cell(length(NDs), length(durs));
    fictracFieldStrct.repsOptoTimes = cell(length(NDs), length(durs));
    fictracFieldStrct.repsPDataNames = cell(length(NDs), length(durs));
    fictracFieldStrct.means = cell(length(NDs), length(durs));
    fictracFieldStrct.stdErrs = cell(length(NDs), length(durs));

    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        fictracOpto.(thisVarName) = fictracFieldStrct;
    end

    % preallocate struct for leg data
    % for each variable, struct for each rep, mean, and std error
    % reps as #NDs x #durs x #tracked points cell array, with each element
    %  being a matrix of reps x rep length
    % repsOptoTimes as #NDs x #durs x #tracked points cell array, with each
    %  element being a vector of # reps length, tracking opto stim start
    %  time corresponding to each rep
    % repsPDataNames as #NDs x #durs x #tracked points cell array, with
    %  each element being a cell array of # reps length, tracking the name
    %  of the pData file the rep came from
    % means and stdErrs as same size cell array, each element as vector
    legFieldStrct.reps = cell(length(NDs),length(durs),NUM_LEG_PTS);
    legFieldStrct.repsOptoTimes = cell(length(NDs),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.repsPDataNames = cell(length(NDs),length(durs),...
        NUM_LEG_PTS);
    legFieldStrct.means = cell(length(NDs),length(durs),NUM_LEG_PTS);
    legFieldStrct.stdErrs = cell(length(NDs),length(durs),NUM_LEG_PTS);

    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};
        legTrackOpto.(thisVarName) = legFieldStrct;
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
        
        % check that pData has legTrack, FicTrac and Opto data, otherwise, skip
        if (any(contains(pDataVarsNames, 'legTrack')) && ...
                any(contains(pDataVarsNames, 'opto')) && ...
                any(contains(pDataVarsNames,'fictracProc')))

            % load pData
            load(pDataFullPath, 'legTrack','opto','fictracProc');

            % loop through all FicTrac varialbes
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
                [fictracOpto.(thisVarName).reps, ...
                    fictracOpto.(thisVarName).repsOptoTimes, ...
                    fictracOpto.(thisVarName).repsPDataNames,...
                    ftDurTs] = ...
                    extractTrialsOpto(...
                    fictracOpto.(thisVarName).reps, ...
                    fictracOpto.(thisVarName).repsOptoTimes, ...
                    fictracOpto.(thisVarName).repsPDataNames,...
                    thisVarVal, ...
                    fictracProc.t, opto, NDs, durs, bwStimDur, ...
                    pDataName, norm2StimStart);
            end
            
            % loop through all legTrack variables
            for j = 1:length(legTrackVarNames)
                thisVarName = legTrackVarNames{j};
                thisVarVal = legTrack.(thisVarName);

                % loop through all tracked points, get reps
                for k = 1:NUM_LEG_PTS
                    [legTrackOpto.(thisVarName).reps(:,:,k), ...
                        legTrackOpto.(thisVarName).repsOptoTimes(:,:,k),...
                        legTrackOpto.(thisVarName).repsPDataNames(:,:,k),...
                        legDurTs] = ...
                        extractTrialsOpto(...
                        legTrackOpto.(thisVarName).reps(:,:,k), ...
                        legTrackOpto.(thisVarName).repsOptoTimes(:,:,k),...
                        legTrackOpto.(thisVarName).repsPDataNames(:,:,k),...
                        thisVarVal(:,k), legTrack.t, opto, NDs, durs, ...
                        bwStimDur, pDataName, false);
                end
            end

            % clear this trial's loaded data
            clear fictracProc opto legTrack
        end
    end

    % compute means and std errs for FicTrac vars
    for i = 1:length(fictracVarNames)
        thisVarName = fictracVarNames{i};
        [fictracOpto.(thisVarName).means, ...
            fictracOpto.(thisVarName).stdErrs] = ...
            computeMeansStdErrsFictracOpto(fictracOpto.(thisVarName).reps,false);
    end

    % compute means and std errs for leg tracking vars
    for i = 1:length(legTrackVarNames)
        thisVarName = legTrackVarNames{i};

        % loop through all tracked points
        for j = 1:NUM_LEG_PTS
            [legTrackOpto.(thisVarName).means(:,:,j), ...
                legTrackOpto.(thisVarName).stdErrs(:,:,j)] = ...
                computeMeansStdErrsFictracOpto(...
                legTrackOpto.(thisVarName).reps(:,:,j),false);
        end
    end

    % add timing vectors to fictracOpto and legTrackOpto structs
    fictracOpto.durTs = ftDurTs;
    legTrackOpto.durTs = legDurTs;

    % full path to save data
    saveFullPath = [savePath filesep flyName '_avgLegFictracOpto.mat'];

    % save data
    save(saveFullPath, 'fictracOpto','legTrackOpto','NDs','durs',...
        'bwStimDur', '-v7.3');

    fprintf('Saved legFictracOpto for %s!\n', flyName);
end