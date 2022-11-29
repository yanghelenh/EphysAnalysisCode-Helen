% plotCondBinnedSpikerateFicTrac.m
%
% Function that plots spike rate vs. a FicTrac variable (yaw, forward, or
%  lateral velocity), conditioned on bins of another FicTrac variable
% Pools over pData files, selected through GUI
% Only on moving bouts
% Wrapper for genCondBinnedSpikeRateFictracPlot()
%
% INPUTS:
%   tDelay - time offset b/w spike rate and FicTrac variables, in sec;
%       negative values are ephys before behavior
%   ftVar1Name - x-axis variable name (as 'yaw', 'fwd', or 'lat')
%   ftVar2Name - y-axis variable name (as 'yaw', 'fwd', or 'lat')
%   ftVar1Bins - vector of length 3: [lower lim, upper lim, number bins]
%       for ftVar1
%   ftVar2Bins - matrix of size m x 2, where each row is [lower lim, upper
%       lim] for one bin; number of bins = m; will produce m separate lines
%   yScale - vector of length 2 for yScale of plot [lower upper]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 10/26/22 - HHY
%
% UPDATED: 10/26/22 - HHY
%

function plotCondBinnedSpikerateFicTrac(tDelay, ftVar1Name, ftVar2Name, ...
    ftVar1Bins, ftVar2Bins, yScale)


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
    % initialize vectors for tracking step info across pData files
    allFtVar1 = [];
    allFtVar2= [];
    allSpikeRate = []; 

    % loop through all pData files
    for i = 1:numPDataFiles 
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];

        % get variables in pData file
        pDataMatObj = matfile(pDataFullPath);
        pDataVarsStrct = whos(pDataMatObj);
        pDataVars = struct2cell(pDataVarsStrct);
        pDataVarsNames = pDataVars(1,:); % cell array of names    

        if (any(contains(pDataVarsNames, 'fictracProc')) &&...
                any(contains(pDataVarsNames, 'ephysSpikes')) && ...
                any(contains(pDataVarsNames, 'moveNotMove')))

            load(pDataFullPath, 'fictracProc', 'ephysSpikes', ...
                'moveNotMove');

            % get ftVar1 values
            if strcmpi(ftVar1Name, 'yaw')
                ftVar1 = fictracProc.yawAngVel(moveNotMove.ftMoveInd);
                xAxisName = 'Yaw velocity (deg/s)';
                ftVar1FullName = 'Yaw Velocity';
            elseif strcmpi(ftVar1Name, 'fwd')
                ftVar1 = fictracProc.fwdVel(moveNotMove.ftMoveInd);
                xAxisName = 'Forward velocity (mm/s)';
                ftVar1FullName = 'Forward Velocity';
            elseif strcmpi(ftVar1Name, 'lat')
                ftVar1 = fictracProc.slideVel(moveNotMove.ftMoveInd);
                xAxisName = 'Lateral velocity (mm/s)';
                ftVar1FullName = 'Lateral Velocity';
            end

            % make sure ftVar1 is column vector
            if ~iscolumn(ftVar1)
                ftVar1 = ftVar1';
            end

            % get ftVar2 values
            if strcmpi(ftVar2Name, 'yaw')
                ftVar2 = fictracProc.yawAngVel(moveNotMove.ftMoveInd);
                ftVar2Units = 'deg/s';
                ftVar2FullName = 'Yaw Velocity';
            elseif strcmpi(ftVar2Name, 'fwd')
                ftVar2 = fictracProc.fwdVel(moveNotMove.ftMoveInd);
                ftVar2Units = 'mm/s';
                ftVar2FullName = 'Yaw Velocity';
            elseif strcmpi(ftVar2Name, 'lat')
                ftVar2 = fictracProc.slideVel(moveNotMove.ftMoveInd);
                ftVar2Units = 'mm/s';
                ftVar2FullName = 'Lateral Velocity';
            end

            % make sure ftVar2 is column vector
            if ~iscolumn(ftVar2)
                ftVar2 = ftVar2';
            end

            % get spikeRate, introduce tDelay
            spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
                fictracProc.t);
            spikeRate = spikeRate(moveNotMove.ftMoveInd);

            % make sure spikeRate is column vector
            if ~iscolumn(spikeRate)
                spikeRate = spikeRate';
            end

            % add to running vectors
            allFtVar1 = [allFtVar1; ftVar1];
            allFtVar2 = [allFtVar2; ftVar2];
            allSpikeRate = [allSpikeRate; spikeRate];
        end
    end

    % remove all NaN spike rates
    notNanLog = ~isnan(allSpikeRate);
    allFtVar1 = allFtVar1(notNanLog);
    allFtVar2 = allFtVar2(notNanLog);
    allSpikeRate = allSpikeRate(notNanLog);

    % feed into plotting function
    [f, meanSpikeRate, stdDevSpikeRate, semSpikeRate, binMids] = ...
    genCondBinnedSpikerateFictracPlot(allSpikeRate, allFtVar1, allFtVar2, ...
        ftVar1Bins, ftVar2Bins, yScale, xAxisName, ...
        'Spike Rate (spikes/s)', ftVar1FullName, ftVar2FullName, ftVar2Units);

end