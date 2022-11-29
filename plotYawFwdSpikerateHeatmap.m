   
function plotYawFwdSpikerateHeatmap(tDelay, spikeRateScale)
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
    allYaw = [];
    allFwd = [];
    allSpikeRate = []; 

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        flyName = pDataName(1:(end-12));
        
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

            fwdVel = fictracProc.fwdVel(moveNotMove.ftMoveInd);
            yawVel = fictracProc.yawAngVel(moveNotMove.ftMoveInd);

            spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
                fictracProc.t);
            spikeRate = spikeRate(moveNotMove.ftMoveInd);

            % add to running vectors
            allYaw = [allYaw; yawVel];
            allFwd = [allFwd; fwdVel];
            allSpikeRate = [allSpikeRate; spikeRate];

        end
    end

    xDataName = 'yawAngVel';
    yDataName = 'fwdVel';
    zDataName = 'spikeRate';

    [f, heatmapMat, countsMat] = genHeatmap(allYaw, allFwd, allSpikeRate,...
        xDataName, yDataName, zDataName, [-500 500 30], [-5 15 30], spikeRateScale, 20,...
        0, [], flyName);
end