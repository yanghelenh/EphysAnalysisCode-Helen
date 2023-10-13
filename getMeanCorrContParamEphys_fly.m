function allFliesCorr = getMeanCorrContParamEphys_fly(datDir, xRange, isCirc)

% prompt user to select getCorrelationEphysContParam() files
        [corrEphysFNames, corrEphysPath] = uigetfile('*.mat', ...
            'Select output file', datDir, 'MultiSelect', 'on');   


    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(corrEphysFNames))
        numFlies = length(corrEphysFNames);
    else
        numFlies = 1;
    end

        % preallocate
    allFliesCorr = zeros(numFlies,1);


    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(corrEphysFNames))
            outName = corrEphysFNames{i};
        else
            outName = corrEphysFNames;
        end

        % load data from file 1
        fullFilePath = [corrEphysPath filesep outName];
    
        % load variables
        load(fullFilePath, 'ephysValsNorm', 'behVals1D', ...
            'ephysParam', 'behParams', 'legs', 'tDelay');    


        % get logical for which samples fall into valid xRange
        valLog = (ephysValsNorm >= xRange(1)) & ...
            (ephysValsNorm < xRange(2));
        % get ephys values for these samples
        thisEphysVal = ephysValsNorm(valLog);
        thisBehVal = behVals1D(valLog);

        % if circular behavioral parameter
        if isCirc
            [thisR, ~] = circ_corrcl(deg2rad(thisBehVal), thisEphysVal);
        else
            % get Pierson correlation coefficient
            thisR = corr(thisEphysVal, thisBehVal);
%             if contains(outName,'220907')
%                 thisR = thisR * -1;
%             end
        end

        allFliesCorr(i) = thisR;

    end
    meanCorr = mean(allFliesCorr);
    semCorr = std(allFliesCorr) / sqrt(numFlies);

    fprintf('r = %0.2f +/- %0.2f\n', meanCorr,semCorr);

end