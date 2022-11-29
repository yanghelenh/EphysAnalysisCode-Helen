% plotAEPPEPbyMultiFictrac.m
%
% Function to plot mean and standard deviation of leg AEPs or PEPs (user 
%  specified), colored by values of two FicTrac variables (yaw, forward, or 
%  lateral velocity; user specified).
% User specifies number of bins in which to divide the FicTrac variables or 
%  the edges of the bins. If number of bins specified, divide spike rate 
%  into that number of bins by percentiles (i.e. equal size bins by number 
%  of steps). Specify bin edges as n x 2 matrix, where n = num bins, and 
%  the first column is the bin start and the second the bin end.
% Works on multiple pData files. Select through GUI. Should be for 1 fly.
% For pData files with iInj or opto, only uses steps during no stim period.
%
% INPUTS:
%   whichEP - which extreme position to plot (AEP or PEP), specified as
%     string 'AEP' or 'PEP'
%   whichFicTrac - which FicTrac variables, specified as cell array of 
%       strings 'yaw', 'fwd', or 'lat'
%   bins - as cell array equal to number of FicTrac variables. For each, 
%       if scalar, number of bins. If n x 2 matrix, n bins with first
%       column as bin starts and second as bin ends
%   whichPhase - which phase of leg movement, specified as string 'swing'
%     or 'stance'
%   
% OUTPUTS:
%   none - but produces plot
%
% CREATED: 10/26/22 - HHY
%
% UPDATED:
%   10/26/22
%
function plotAEPPEPbyMultiFictrac(whichEP, whichFicTrac, bins, whichPhase)
    
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
    allStepFictrac = cell(size(whichFicTrac)); % FicTrac val for each step

    % number of FicTrac variables to bin over
    numVars = length(whichFicTrac);

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
        if (any(contains(pDataVarsNames, 'legSteps')))

            % load pData: if current injection, load that info
            % also load appropriate set of step parameters based on
            %  swing/stance phase
            switch lower(whichPhase)
                case 'swing'
                    if (any(contains(pDataVarsNames, 'legStepsByIinj')))
                        load(pDataFullPath, 'legSteps', 'legStepsByIinj',...
                            'swingStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByIinj.stepIinjCat, ...
                            legSteps.stepSwingStance,-1);
                    elseif (any(contains(pDataVarsNames, 'legStepsByOpto')))
                        load(pDataFullPath, 'legSteps', 'legStepsByOpto',...
                            'swingStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByOpto.stepOptoCat, ...
                            legSteps.stepSwingStance,-1);
                    else
                        load(pDataFullPath, 'swingStepParams');
                    end
                    % name step parameters for this phase name
                    stepParams = swingStepParams;
                case 'stance'
                    if (any(contains(pDataVarsNames, 'legStepsByIinj')))
                        load(pDataFullPath, 'legSteps', 'legStepsByIinj',...
                            'stanceStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByIinj.stepIinjCat, ...
                            legSteps.stepSwingStance,1);

                    elseif (any(contains(pDataVarsNames, 'legStepsByOpto')))
                        load(pDataFullPath, 'legSteps', 'legStepsByOpto',...
                            'stanceStepParams');

                        % get step categories for this phase
                        thisLegStepCat = groupStepParamsBySwingStance(...
                            legStepsByOpto.stepOptoCat, ...
                            legSteps.stepSwingStance,1);
                    else
                        load(pDataFullPath,'stanceStepParams');
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
            % get step leg info
            stepWhichLeg = stepParams.stepWhichLeg;

            % get step FicTrac info
            stepFt = cell(size(whichFicTrac));
            ftUnits = cell(size(whichFicTrac));
            ftName = cell(size(whichFicTrac));

            for j = 1:numVars
                switch lower(whichFicTrac{j})
                    case 'yaw'
                        stepFt{j} = stepParams.stepFtYaw;
                        ftUnits{j} = 'deg/s'; % units for plotting
                        ftName{j} = 'FicTrac Yaw Velocity'; % for plotting
                    case 'fwd'
                        stepFt{j} = stepParams.stepFtFwd;
                        ftUnits{j} = 'mm/s'; % units for plotting
                        ftName{j} = 'FicTrac Forward Velocity'; % for plotting
                    case 'lat'
                        stepFt{j} = stepParams.stepFtLat;
                        ftUnits{j} = 'mm/s'; % units for plotting
                        ftName{j} = 'FicTrac Lateral Velocity'; % for plotting
                end
            end

            % if this pData file has current injection, extract only steps
            %  in no stim period
            if (any(contains(pDataVarsNames, 'legStepsByIinj')) || ...
                    any(contains(pDataVarsNames, 'legStepsByOpto')))
                stepEPX = stepEPX(thisLegStepCat == NO_STIM_IND);
                stepEPY = stepEPY(thisLegStepCat == NO_STIM_IND);

                stepWhichLeg = stepWhichLeg(thisLegStepCat == NO_STIM_IND);

                for j = 1:length(whichFicTrac)
                    stepFt{j} = stepFt{j}(thisLegStepCat == NO_STIM_IND);
                end
            end

            % append these steps to running vector across pData files
            allStepXvals = [allStepXvals; stepEPX];
            allStepYvals = [allStepYvals; stepEPY];
            allStepLegs = [allStepLegs; stepWhichLeg];
            for j = 1:length(whichFicTrac)
                allStepFictrac{j} = [allStepFictrac{j}; stepFt{j}];
            end
        end
    end

    % get bins for grouping EPs by the FicTrac variables
    numBins = 1; % start as default, all data together
    binStarts = cell(size(whichFicTrac));
    binEnds = cell(size(whichFicTrac));
    % number of bins for each FicTrac variable
    numBinsEach = zeros(size(whichFicTrac));

    for j = 1:numVars
        if (length(bins{j}) == 1) % if scalar value for number of bins
            % number of bins for this FicTrac var
            thisNumBins = bins{j};
            numBinsEach(j) = thisNumBins;
    
            % total number of bins is multiplied across all FicTrac vars
            numBins = numBins * thisNumBins;
            
            % get bin edges, using percentiles
            prcts = linspace(0,100,thisNumBins + 1); % bin edges, as percentiles
            % convert percentiles to units of the FicTrac variable
            binEdges = prctile(allStepFictrac{j},prcts);
            % convert edges to starts and ends, as separate vectors
            binStarts{j} = binEdges(1:(end-1));
            binEnds{j} = binEdges(2:end);
    
        else % bin edges specified in function input
            thisNumBins = size(bins{j},1);
            numBinsEach(j) = thisNumBins;
    
            numBins = numBins * thisNumBins;
            binStarts{j} = bins{j}(:,1);
            binEnds{j} = bins{j}(:,2);
        end
    end

    % map number of bins to every bin combination across FicTrac variables
    % works for 1-3 FicTrac variables
    binInds = cell(size(whichFicTrac));
    for j = 1:numVars
        switch j
            case 1
                binInds{j} = repmat(1:numBinsEach(j),1,numBins / ...
                    numBinsEach(j));
            case 2
                binInds{j} = repelem(1:numBinsEach(j),numBins / ...
                    numBinsEach(j));
            case 3
                binInds{j} = repmat(repelem(1:numBinsEach(j),numBinsEach(j)),...
                    1, numBins/(numBinsEach(j)^2));
        end
    end

    % preallocate, matrices for means and stddev
    xMeans = zeros(numBins, NUM_LEGS);
    yMeans = zeros(numBins, NUM_LEGS);
    xStdDev = zeros(numBins, NUM_LEGS);
    yStdDev = zeros(numBins, NUM_LEGS);

    % get means and std dev, by leg
    for i = 1:numBins
        % get indices of steps that belong in this bin
        switch numVars
            case 1
                thisBinInds = find(...
                    (allStepFictrac{1} >= binStarts{1}(binInds{1}(i))) & ...
                    (allStepFictrac{1} < binEnds{1}(binInds{1}(i))) ...
                    );
            case 2
                thisBinInds = find(...
                    (allStepFictrac{1} >= binStarts{1}(binInds{1}(i))) & ...
                    (allStepFictrac{1} < binEnds{1}(binInds{1}(i))) & ...
                    (allStepFictrac{2} >= binStarts{2}(binInds{2}(i))) & ...
                    (allStepFictrac{2} < binEnds{2}(binInds{2}(i))) ...
                    );
            case 3
                thisBinInds = find(...
                    (allStepFictrac{1} >= binStarts{1}(binInds{1}(i))) & ...
                    (allStepFictrac{1} < binEnds{1}(binInds{1}(i))) & ...
                    (allStepFictrac{2} >= binStarts{2}(binInds{2}(i))) & ...
                    (allStepFictrac{2} < binEnds{2}(binInds{2}(i))) & ...
                    (allStepFictrac{3} >= binStarts{3}(binInds{3}(i))) & ...
                    (allStepFictrac{3} < binEnds{3}(binInds{3}(i)))...
                    );
        end
        
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

        switch numVars
            case 1    
                legendStr{i} = sprintf('%s: %.f to %.f %s', ...
                    whichFicTrac{1}, binStarts{1}(binInds{1}(i)), ...
                    binEnds{1}(binInds{1}(i)), ftUnits{1});
            case 2
                legendStr{i} = sprintf('%s: %.f to %.f %s, %s: %.f to %.f %s', ...
                    whichFicTrac{1}, binStarts{1}(binInds{1}(i)), ...
                    binEnds{1}(binInds{1}(i)), ftUnits{1}, ...
                    whichFicTrac{2}, binStarts{2}(binInds{2}(i)), ...
                    binEnds{2}(binInds{2}(i)), ftUnits{2} ...
                    );
            case 3
                legendStr{i} = sprintf(...
                    '%s: %.f to %.f %s, %s: %.f to %.f %s, %s: %.f to %.f %s', ...
                    whichFicTrac{1}, binStarts{1}(binInds{1}(i)), ...
                    binEnds{1}(binInds{1}(i)), ftUnits{1}, ...
                    whichFicTrac{2}, binStarts{2}(binInds{2}(i)), ...
                    binEnds{2}(binInds{2}(i)), ftUnits{2}, ...
                    whichFicTrac{3}, binStarts{3}(binInds{3}(i)), ...
                    binEnds{3}(binInds{3}(i)), ftUnits{3} ...
                    );
        end
    end
    % x and y lims
    xlim([-1 1]);
    ylim([-1 1]);

    % x and y labels
    xlabel('Body Lengths <-L - R->');
    ylabel('Body Lengths')

    % title of plot: EP and time delay in ms
    ttlStr = sprintf('%s', whichEP);
    title(ttlStr);

    % reverse y axis (x values) so head (neg vals) is at top
    set(gca, 'YDir','reverse');

    legend(legendStr); % legend 
end