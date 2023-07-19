% computeAvgFictracOpto_1Fly.m
%
% Function that extracts the mean, std err, and individual responses
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
%   flipLR - boolean for whether to flip left-right (true when stim cell is
%       on left)
%   savePath - path to folder in which to save data
%
% OUTPUTS:
%   none, but saves struct to file of name date_fly_avgFictracOpto.mat into
%       savePath
%
% CREATED: 3/29/22 - HHY
%
% UPDATED:
%   3/29/22 - HHY
%   7/10/23 - HHY - updated with flipping left/right
%
function computeAvgFictracOpto_1Fly(NDs, durs, bwStimDur, flipLR, savePath)

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

    % preallocate
    fwdVelReps = cell(length(NDs), length(durs));
    slideVelReps = cell(length(NDs), length(durs));
    yawVelReps = cell(length(NDs), length(durs));
    
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
        
        % load pData
        load(pDataFullPath, 'exptCond');
        
        % check that pData has FicTrac and Opto data, otherwise, skip
        if ((contains(exptCond, 'Fictrac', 'IgnoreCase', true)) && ...
                (contains(exptCond, 'Opto', 'IgnoreCase', true)))

            % load pData
            load(pDataFullPath, 'fictracProc','opto')

            % for each FicTrac variable, replace dropped indices with NaN
            fwdVel = fictracProc.fwdVel;
            fwdVel(fictracProc.dropInd) = nan;
            slideVel = fictracProc.slideVel;
            slideVel(fictracProc.dropInd) = nan;
            yawVel = fictracProc.yawAngVel;
            yawVel(fictracProc.dropInd) = nan;

            % flip yaw and slide if flipping left/right
            if (flipLR)
                slideVel = slideVel * -1;
                yawVel = yawVel * -1;
            end

            % get all reps for this pData for each FicTrac variable
            [fwdVelReps,~] = extractFictracTrialsOpto(fwdVelReps, fwdVel, ...
                fictracProc.t, opto, NDs, durs, bwStimDur);
            [slideVelReps,~] = extractFictracTrialsOpto(slideVelReps, slideVel, ...
                fictracProc.t, opto, NDs, durs, bwStimDur);
            [yawVelReps, durTs] = extractFictracTrialsOpto(yawVelReps, yawVel, ...
                fictracProc.t, opto, NDs, durs, bwStimDur);

            % clear this trial's loaded data
            clear fictracProc opto exptCond
        end
    end

    % compute means and std errs
    [fwdVelMeans, fwdVelStdErrs] = ...
        computeMeansStdErrsFictracOpto(fwdVelReps);
    [slideVelMeans, slideVelStdErrs] = ...
        computeMeansStdErrsFictracOpto(slideVelReps);
    [yawVelMeans, yawVelStdErrs] = ...
        computeMeansStdErrsFictracOpto(yawVelReps);

    % put data in struct
    fwdVelData.means = fwdVelMeans;
    fwdVelData.stdErrs = fwdVelStdErrs;
    fwdVelData.reps = fwdVelReps;
    slideVelData.means = slideVelMeans;
    slideVelData.stdErrs = slideVelStdErrs;
    slideVelData.reps = slideVelReps;
    yawVelData.means = yawVelMeans;
    yawVelData.stdErrs = yawVelStdErrs;
    yawVelData.reps = yawVelReps;

    % full path to save data
    saveFullPath = [savePath filesep flyName '_avgFictracOpto.mat'];

    % save data
    save(saveFullPath, 'fwdVelData','slideVelData','yawVelData', 'NDs',...
        'durs', 'durTs', 'flipLR', '-v7.3');

end