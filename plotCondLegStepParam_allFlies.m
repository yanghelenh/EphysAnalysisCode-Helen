% plotCondLegStepParam_allFlies.m
%
% Function that takes output files from saveLegStepParamByCond_fly() and
%  generates a plot of the legStep param for each leg. Plots mean +/- SEM 
%  for each fly. Different conditions in different columns. Points from the
%  same fly are connected by lines. 
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   yScale - scale for plots, as [min max]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 2/14/23 - HHY
%
% UPDATED:
%   2/14/23 - HHY
%
function plotCondLegStepParam_allFlies(datDir, yScale)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from saveAEPPEPbyCond_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select Step Param files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % initialize figure
    figure;
%     c = colormap('lines');

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'stepValMeans', 'stepValSEM', 'legStepParam',...
            'stim', 'flipLegsLR');

        if (flipLegsLR)
            if (strcmpi(legStepParam, 'stepDirections'))
                stepValMeans = -1 * stepValMeans;
            end
        end

        % number of conditions
        numCats = size(stepValMeans,1);
        
        % x vector for plotting (for number of conditions)
        xVec = 1:numCats;

        for j = 1:6
            subplot(3,2,subInd(j));
            hold on;

            % plot for this fly
            errorbar(xVec', stepValMeans(:,j), stepValSEM(:,j), ...
                'Marker', 'x','LineWidth',1, 'CapSize', 0, 'Color', 'k');
    
            ylim(yScale);
            xticks(xVec);
            xticklabels(repmat({''},1,numCats));
        end
    end

    

    % generate labels for x axis
    xLabels = cell(size(xVec));
    xLabels{1} = 'No Stim';

    if (strcmpi(stim.whichStim, 'iInj'))
        % generate key for mapping b/w indices and amps and durs
        iInjCatAmps = zeros(1,numCats);
        iInjCatDurs = zeros(1,numCats);
        % counter index into vectors, skip 1 for 0,0, no iInj
        counter = 2; 
    
        % assign amps to indices
        for i = 1:length(stim.amps)
            for j = 1:length(stim.durs)
                iInjCatAmps(counter) = stim.amps(i);
                iInjCatDurs(counter) = stim.durs(j);
    
                counter = counter + 1;
            end
        end
    elseif (strcmpi(stim.whichStim, 'opto'))
        % generate key for mapping b/w indices and NDs and durs
        optoCatNDs = ones(1,numCats) * -1;
        optoCatDurs = ones(1,numCats) * -1;
        % counter index into vectors, skip 1 for -1, -1, no opto
        counter = 2; 

        % assign NDs, durs to indices
        for i = 1:length(stim.NDs)
            for j = 1:length(stim.durs)
                optoCatNDs(counter) = stim.NDs(i);
                optoCatDurs(counter) = stim.durs(j);

                counter = counter + 1;
            end
        end
    end

    for i = 2:numCats
        if (strcmpi(stim.whichStim, 'iInj'))
            xLabels{i} = sprintf('%d pA, %.1f s', iInjCatAmps(i), ...
                iInjCatDurs(i));
        elseif (strcmpi(stim.whichStim, 'opto'))
            xLabels{i} = sprintf('ND=%.1f, %.1f s', optoCatNDs(i), ...
                optoCatDurs(i));
        end
    end

    xticklabels(xLabels);


    sgtitle(legStepParam);
end
