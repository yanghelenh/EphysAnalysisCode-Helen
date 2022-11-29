% plotAEPPEPbySpikeRate.m
%
% Function to plot mean and standard deviation of leg AEPs or PEPs (user 
%  specified), colored by spike rate.
% User specified delay between ephys and behavior. Negative times are ephys
%  before behavior.
% User specifies number of bins in which to divide spike rate or the edges
%  of the bins. If number of bins specified, divide spike rate into that
%  number of bins by percentiles (i.e. equal size bins by number of steps).
%  Specify bin edges as n x 2 matrix, where n = num bins, and the first
%  column is the bin start and the second the bin end.
% Works on multiple pData files. Select through GUI. Should be for 1 fly.
% For pData files with iInj, only uses steps during no stim period.
%
% INPUTS:
%   whichEP - which extreme position to plot (AEP or PEP), specified as
%     string 'AEP' or 'PEP'
%   bins - if scalar, number of bins. If n x 2 matrix, n bins with first
%       column as bin starts and second as bin ends
%   tDelay - time (in sec) to shift spike and step times relative to each
%       other. Negative times are ephys before behavior. Scalar
%   whichPhase - which phase of leg movement, specified as string 'swing'
%     or 'stance'
%   
% OUTPUTS:
%   none - but produces plot
%
% CREATED: 10/7/22 - HHY
%
% UPDATED:
%   10/7/22
%
function plotAEPPEPbySpikeRate(whichEP, bins, tDelay, whichPhase)
    
    legInd = 1:6; % indicies into raw position matricies for legs
    NUM_LEGS = 6; % number of legs

    NO_STIM_IND = 1; % code for iInj corresponding to no stim

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
    allStepXvals = [];
    allStepYvals = [];
    allStepLegs = []; % which leg each step belongs to
    allStepSpikeRate = []; % spike rate for each step

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

        % check that pData has appropriate structs; otherwise, skip
        if (any(contains(pDataVarsNames, 'legSteps')) && ...
                any(contains(pDataVarsNames, 'ephysSpikes')))

            % load pData: if current injection, load that info
            % also load appropriate set of step parameters based on
            %  swing/stance phase
            switch lower(whichPhase)
                case 'swing'
                    if (any(contains(pDataVarsNames, 'legStepsByIinj')))
                        load(pDataFullPath, 'legSteps', 'legStepsByIinj',...
                            'ephysSpikes', 'swingStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByIinj.stepIinjCat, ...
                            legSteps.stepSwingStance,-1);
                    else
                        load(pDataFullPath, 'ephysSpikes', ...
                            'swingStepParams');
                    end
                    % name step parameters for this phase name
                    stepParams = swingStepParams;
                case 'stance'
                    if (any(contains(pDataVarsNames, 'legStepsByIinj')))
                        load(pDataFullPath, 'legSteps', 'legStepsByIinj',...
                            'ephysSpikes', 'stanceStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByIinj.stepIinjCat, ...
                            legSteps.stepSwingStance,1);
                    else
                        load(pDataFullPath, 'ephysSpikes', ...
                            'stanceStepParams');
                    end
                    stepParams = stanceStepParams;
            end

            % get step EP info, for specified EP
            switch lower(whichEP)
                case 'aep'
                    stepEPX = stepParams.stepAEPX;
                    stepEPY = stepParams.stepAEPY;
                case 'pep'
                    stepEPX = stepParams.stepPEPX;
                    stepEPY = stepParams.stepPEPY;
            end
            % get step leg and timing info
            stepWhichLeg = stepParams.stepWhichLeg;
            stepTimes = stepParams.stepT;

            % if this pData file has current injection, extract only steps
            %  in no stim period
            if (any(contains(pDataVarsNames, 'legStepsByIinj')))
                stepEPX = stepEPX(thisLegStepCat == NO_STIM_IND);
                stepEPY = stepEPY(thisLegStepCat == NO_STIM_IND);

                stepWhichLeg = stepWhichLeg(thisLegStepCat == NO_STIM_IND);
                stepTimes = stepTimes(thisLegStepCat == NO_STIM_IND, :);
            end

            % get spike rate for each step
            spikeTimes = ephysSpikes.t(ephysSpikes.startInd);
            stepSpikeRate = getSpikeRateSteps(spikeTimes, stepTimes, ...
                tDelay);

            % append these steps to running vector across pData files
            allStepXvals = [allStepXvals; stepEPX];
            allStepYvals = [allStepYvals; stepEPY];
            allStepLegs = [allStepLegs; stepWhichLeg];
            allStepSpikeRate = [allStepSpikeRate; stepSpikeRate]; 
        end
    end

    % get bins for grouping EPs by spike rate
    if (length(bins) == 1) % if scalar value for number of bins
        numBins = bins;
        
        % get bin edges, using percentiles
        prcts = linspace(0,100,numBins + 1); % bin edges, as percentiles
        % convert percentiles to spike rates
        binEdges = prctile(allStepSpikeRate,prcts);
        % convert edges to starts and ends, as separate vectors
        binStarts = binEdges(1:(end-1));
        binEnds = binEdges(2:end);

    else % bin edges specified in function input
        numBins = size(bins,1);
        binStarts = bins(:,1);
        binEnds = bins(:,2);
    end

    % preallocate, matrices for means and stddev
    xMeans = zeros(numBins, NUM_LEGS);
    yMeans = zeros(numBins, NUM_LEGS);
    xStdDev = zeros(numBins, NUM_LEGS);
    yStdDev = zeros(numBins, NUM_LEGS);

    % get means and std dev, by leg
    for i = 1:numBins
        % get indices of steps that belong in this bin
        thisBinInds = find((allStepSpikeRate >= binStarts(i)) & ...
            (allStepSpikeRate < binEnds(i)));
        % get X, Y, whichLeg for steps in this bin
        thisBinX = allStepXvals(thisBinInds);
        thisBinY = allStepYvals(thisBinInds);
        thisBinLegs = allStepLegs(thisBinInds);

        % loop over all legs, sort steps by leg, then get mean and std dev
        for j = 1:NUM_LEGS
            thisLegX = thisBinX(thisBinLegs == legInd(j));
            xMeans(i,j) = mean(thisLegX);
            xStdDev(i,j) = std(thisLegX);

            thisLegY = thisBinY(thisBinLegs == legInd(j));
            yMeans(i,j) = mean(thisLegY);
            yStdDev(i,j) = std(thisLegY);
        end
    end


    % PLOTTING

    % initialize cell array for legend
    legendStr = cell(1,numBins);

    figure;

    % use colormap lines
    c = colormap('lines');

    for i = 1:numBins
        errorbar(yMeans(i,:), xMeans(i,:), xStdDev(i,:), xStdDev(i,:), ...
            yStdDev(i,:), yStdDev(i,:), ...
            'Marker','x', 'LineStyle','none','Color',c(i,:), ...
            'LineWidth',1.5);

        hold on;

        legendStr{i} = sprintf('%.f Hz to %.f Hz', binStarts(i), ...
            binEnds(i));
    end
    % x and y lims
    xlim([-1 1]);
    ylim([-1 1]);

    % x and y labels
    xlabel('Body Lengths <-L - R->');
    ylabel('Body Lengths')

    % title of plot: EP and time delay in ms
    ttlStr = sprintf('%s binned by spike rate, tDelay = %d ms', whichEP,...
        tDelay * 1000);
    title(ttlStr);

    % reverse y axis (x values) so head (neg vals) is at top
    set(gca, 'YDir','reverse');

    legend(legendStr); % legend 
end