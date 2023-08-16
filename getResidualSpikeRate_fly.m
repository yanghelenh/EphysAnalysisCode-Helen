% getResidualSpikeRate_fly.m
%
% Function for calculating residual spike rate. Takes in input from
%  getCorrelationEphysContParam_cond_fly() and generates residual spike
%  rate. Uses behVals and ephysVals. Takes 1st column of behVals as
%  variable to get mean curve and subtract from this to get residuals.
%  Remaining columns of behVals combined into 1 equally weighted variable
%  and correlated against residuals.
%
% Select getCorrelationEphysContParam_cond_fly() output file through GUI
%
% INPUTS
%   datDir - full path to folder with 
%       getCorrelationEphysContParam_cond_fly() output files
%   behVarRange - 2 element vector with min and max values for explanatory
%       behavioral variable
%   numBins - number of bins to bin explanatory behavioral variable
%   minNumVals - minimum number of values per bin to be included
%   normEphys - boolean for whether to use normalized ephys
%
% OUTPUTS:
%   none, but generates file
%
% CREATED: 8/14/23 - HHY
%
% UPDATED:
%   8/14/23 - HHY
%
function getResidualSpikeRate_fly(datDir, behVarRange, numBins, ...
    minNumVals, normEphys)

    % get input file (output of getCorrelationEphysContParam_cond_fly())
    [inFName, inFDir] = uigetfile('*.mat', ...
        'Select corrEphysParam file', datDir, 'MultiSelect', 'off');

    % load data from this file
    inFileFullPath = [inFDir filesep inFName];

    load(inFileFullPath, 'behVals', 'ephysVals', 'ephysValsNorm');

    % get boundaries of bins
    binSize = (behVarRange(2) - behVarRange(1)) / numBins;
    binEdges = behVarRange(1):binSize:behVarRange(2);
    binStarts = binEdges(1:(end-1));
    binEnds = binEdges(2:end);
    binMids = (binStarts + binEnds)/2;

    % preallocate 
    meanEphysVal = zeros(length(binMids),1);

    % loop through all bins, find mean for each bin
    for j = 1:numBins
        % get logical for which samples fall into this bin
        % first column is variable to bin on
        thisBinLog = (behVals(:,1) >= binStarts(j)) & ...
            (behVals(:,1) < binEnds(j));
        % get ephys values for these samples
        if normEphys
            thisEphys = ephysValsNorm(thisBinLog);
        else
            thisEphys = ephysVals(thisBinLog);
        end

        % check if there are enough time points for this bin
        if (length(thisEphys) >= minNumVals)         
            % get mean and SEM for this bin
            meanEphysVal(j) = mean(thisEphys);
        else
            meanEphysVal(j) = nan;
        end
    end

    % remove NaNs
    nanLog = isnan(meanEphysVal);
    meanEphysVal(nanLog) = [];
    binMids(nanLog) = [];

    % fit curve to bin means
    fitobj = fit(binMids', meanEphysVal, 'exp1');

    % use curve to get expected spike rate for behVals
    expSpikeRate = feval(fitobj, behVals(:,1));

    % get residuals of spike rate
    if normEphys
        residSpikeRate = ephysValsNorm - expSpikeRate;
    else
        residSpikeRate = ephysVals - expSpikeRate;
    end

    % if multiple behavioral variables to explain residual, get 1D equal
    %  weighting
    if (size(behVals,2) > 2)
        % get coefficients
        coeffs = ones(size(behVals,2)-1, 1) * (1/(size(behVals,2)-1));
        residBehVars1D = getLinProj(behVals(:,2:end), coeffs);
    else
        residBehVars1D = behVals(:,2);
    end

    % save output to same file as input, just add additional variables
    save(inFileFullPath, 'residSpikeRate', 'residBehVars1D', ...
        'fitobj', 'expSpikeRate', 'normEphys', '-append');
end