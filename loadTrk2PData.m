% loadTrk2PData.m
%
% Function that finds the .trk file corresponding to the selected pData
%  file(s), preprocesses the .trk data, and saves it into the pData file
% Assumes that .trk file name matches pData name in date_fly#_cell#_trial#
% .trk file as standard output from APT
%
% INPUTS:
%   pDataPath - full path to folder containing pData files
%   trkFilePath - full path to all .trk files
%   trkSuffix - "v#_mdn.trk" - version number and model type and .trk
%       suffix for .trk files; specify in case there are multiple .trk
%       files for a given trial
%   refPts - struct of parameters defining reference tracked points
%       (optional, pass [] to use defaults)
%   smoParams - struct of parameters for smoothing leg postion and 
%       velocity (optional, pass [] to use defaults)
%   
% OUTPUTS:
%   none, but saves processed .trk data into pData file
%
% CREATED: 5/10/22
%
% UPDATED:
%   5/10/22
%
function loadTrk2PData(pDataPath, trkFilePath, trkSuffix, refPts, ...
    smoParams)

    % SOME DEFAULT CONSTANTS %
    % if refPts not specified, default params
    if (isempty(refPts))
        % indicies for specific body parts
        % fit line to head and 3 pts defining midpt b/w legs
        refPts.lineFitInd = 7:10; 
        refPts.midPtInd = 9;
        refPts.headPtInd = 7;
        refPts.abdPtInd = 11;
    end

    % if smoParams not specified, default params
    if (isempty(smoParams))
        % parameters for smoothing for determining leg velocities
        smoParams.padLen = 50; % pad length in samples
        smoParams.sigma = 10; % in samples
       % parameters for removing outliers
       smoParams.medFiltDeg = 5; % degree of median filter
       % percent of (max-min) is acceptable as deviation
       smoParams.percntDev = 10; 
       smoParams.maxMinWin = 1000; % size of window to compute max and min
    end


    % check that trkFilePath does not end with filesep
    if (strcmp(trkFilePath(end),filesep))
        % remove if last character is filesep
        trkFilePath = trkFilePath(1:(end-1));
    end

    % get names of all .trk files in folder specified by trkFilePath with
    %  the trkSuffix
    trkFilesSearch = [trkFilePath filesep '*' trkSuffix];
    trkFilesStrct = dir(trkFilesSearch);

    % convert files struct into cell array of names
    trkFilesNames = cell(size(trkFilesStrct));

    for i = 1:length(trkFilesStrct)
        trkFilesNames{i} = trkFilesStrct(i).name;
    end

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

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        % extract file name without _pData.mat suffix
        % assumes format of date_fly##_cell##_trial##_pData.mat
%         trialFullName = pDataName(1:(end-10));
        trialFullName = pDataName(1:(end-4));

        % full path for this pData file
        pDataFullPath = [pDataDirPath pDataName];
        
        % load only exptCond variable of pData
        load(pDataFullPath, 'exptCond');
        
        % check that pData has leg data, otherwise, skip
        if (contains(exptCond, 'Leg', 'IgnoreCase', true))
            % find .trk file that matches this pData file
            % search by name (date_fly##_cell##_trial##)
            
            % logical of whether trial full name matches names in trk files
            %  list
            trkLog = contains(trkFilesNames, trialFullName);

            % only proceed if there is a matching .trk file
            if (any(trkLog))
                % get .trk file name
                % assumes only 1 matching .trk file (which should be the
                % case)
                matchedTrkFileName = trkFilesNames{trkLog};

                % .trk full path
                trkFileFullPath = [trkFilePath filesep matchedTrkFileName];
                
                % get timing info (leg.frameTimes)
                load(pDataFullPath, 'leg');

                % preprocess .trk data
                legTrack = preprocessLegTrack(trkFileFullPath, ...
                    leg.frameTimes, refPts, smoParams);

                % save legTrack struct into pData
                save(pDataFullPath, 'legTrack', '-append'); 

                % display
                fprintf('Saved legTrack for %s!\n', pDataName);
            else
                fprintf('%s does not have a matching .trk file\n', ...
                    pDataName);
            end
        end
    end
end