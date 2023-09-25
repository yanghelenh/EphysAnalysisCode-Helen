% plotPairCorrContParamEphys_fly.m
%
% Function to plot correlation coefficient between ephys and behavior, 
%  paired by same fly, 2 conditions. Operates on output of 
%  getCorrelationEphysContParam_cond_all().
% Prompts user to select output files 2x, one for each condition. Assumes
%  output files for same fly selected 
%
% INPUTS:
%   datDir - path to folder containing savePCAinterpStepParam() output files
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/25/23 - HHY
%
% UPDATED:
%   9/25/23 - HHY
%
function allFliesCorr = plotPairCorrContParamEphys_fly(datDir, xRange, yScale)

    for n = 1:2
        % prompt user to select getCorrelationEphysContParam() files
        [corrEphysFNames{n}, corrEphysPath{n}] = uigetfile('*.mat', ...
            'Select output file', datDir, 'MultiSelect', 'on');
    
        % sort
        corrEphysFNames{n} = sort(corrEphysFNames{n});
    end

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(corrEphysFNames{1}))
        numFlies = length(corrEphysFNames{1});
    else
        numFlies = 1;
    end

    % preallocate
    allFliesCorr = zeros(numFlies, 2);


    for i = 1:numFlies
        for n = 1:2
            % handle whether it's a cell array or not
            if (iscell(corrEphysFNames{n}))
                outName = corrEphysFNames{n}{i};
            else
                outName = corrEphysFNames{n};
            end
    
            % load data from file 1
            fullFilePath = [corrEphysPath{n} filesep outName];
        
            % load variables
            load(fullFilePath, 'ephysValsNorm', 'behVals1D', ...
                'ephysParam', 'behParams', 'legs', 'tDelay');

    
            % get logical for which samples fall into valid xRange
            valLog = (ephysValsNorm >= xRange(1)) & ...
                (ephysValsNorm < xRange(2));
            % get ephys values for these samples
            thisEphysVal = ephysValsNorm(valLog);
            thisBehVal = behVals1D(valLog);
    
            % get Pierson correlation coefficient
            thisR = corr(thisEphysVal, thisBehVal);
    
            allFliesCorr(i,n) = thisR;

        end
    end


    % x vals for plotting
    xVals = 1:2;
    xScale = [0.5 2.5];

    % generate figure
    figure;

    hold on;

    totAvg = zeros(1,2);
    % get average across flies

    for i = 1:size(allFliesCorr,2)
        thisCond = allFliesCorr(:,i);
        totAvg(i) = mean(thisCond);
    end

    % plot avg across flies
    plot(xVals,totAvg,'Color','k','LineWidth',2, 'Marker', '_');

    % plot individual flies
    plot(xVals,allFliesCorr, 'Marker','.');


    behLblStr = [];
    if iscell(behParams)
        for i = 1:length(behParams)
            if ~i==1
                behLblStr = [behLblStr ', ' behParams{i} ' ' legs{i}];
            else
                behLblStr = [behParams{i} ' ' legs{i}];
            end
        end
    else
        behLblStr = [behParams ' ' legs];
    end
    ylabel(behLblStr);

    xlim(xScale);
    ylim(yScale);
end