% plotLegFictracIndivOptoTrials.m
%
% Function to plot leg and FicTrac data for individual opto stimulation
%  trials for a single fly. Each column is an opto stim trial, each row is
%  a different parameter.
% User specifies which leg and FicTrac parameters, by name of fields in
%  legTrackOpto and fictracOpto structs
% Data expected: output from computeLegFictracOpto_1Fly(). Load in through
%  GUI
%
% INPUTS:
%   paramNames - cell array of names of leg and FicTrac parameters to plot,
%       must match names of fields of legTrackOpto and fictracOpto structs
%   dataDir - directory in which to look for data
%   whichTrials - vector of indices of which trials to plot; invalid
%       indices are ignored
%   whichND - scalar value, which ND to plot
%   whichDur - scalar value, which duration of stimulation to plot
%   legInd - vector of indices of leg points, in order [R1 R2 R3 L1 L2 L3]
%   yScale - y axis scale for all parameters, as n parameters by 2 matrix,
%       where column 1 is min and 2 is max
%
% OUTPUTS:
%   none, but produces plot
%
% CREATED: 5/16/22 - HHY
%
% UPDATED:
%   5/16/22 - HHY
%
function plotLegFictracIndivOptoTrials(paramNames, dataDir, whichTrials,...
    whichND, whichDur, legInd, yScale)

    % smoothing parameter for FicTrac variables, moving average filter
    avgWin = 0.15;

    % prompt user to select avgLegFictracOpto file
    [dataFileName, dataFilePath] = uigetfile('*_avgLegFictracOpto.mat', ...
        'Select avgLegFictracOpto file', dataDir, 'MultiSelect', 'off');

    % load avgLegFictracOpto file
    load([dataFilePath filesep dataFileName], 'legTrackOpto', ...
        'fictracOpto', 'durs', 'NDs', 'bwStimDur');

    % convert whichND and whichDur to indices into legTrackOpto and
    %  fictracOpto cell arrays
    ndInd = find(whichND == NDs);
    durInd = find(whichDur == durs);

    % number of parameters to plot - becomes number of rows in plot
    numParams = length(paramNames);

    % number of trials - becomes number of columns in plot
    numTrials = length(whichTrials);

    % for each parameter, determine whether it's a leg or a FicTrac
    %  parameter
    % initialize cell array to keep track of this info
    paramType = cell(size(paramNames));
    % loop through all parameters
    for i = 1:numParams
        if (isfield(fictracOpto, paramNames{i}))
            paramType{i} = 'fictrac';
        elseif (isfield(legTrackOpto, paramNames{i}))
            paramType{i} = 'legTrack';
        end
    end

    % initialize figure
    figure;
    % set figure size and position, [x y w h] (these numbers for laptop)
    set(gcf,'Position',[10 10 1610 930]);

    % loop through all trials
    for i = 1:numTrials
        % initialize variables to match trials
        thisTrialOptoTime = [];
        thisTrialPDatName = [];

        % initialize vector to hold axes handles
        axHnd = zeros(size(paramNames));

        % loop through all parameters
        for j = 1:numParams
            % which trial is defined by index in first parameter
            if (j == 1)
                firstTrialInd = whichTrials(i);
                thisTrial = firstTrialInd;
                % get opto stim time and pData file name for this trial
                switch paramType{j}
                    case 'fictrac'
                        thisTrialOptoTime = fictracOpto.(paramNames{j}).repsOptoTimes{ndInd,durInd}(firstTrialInd);
                        thisTrialPDatName = fictracOpto.(paramNames{j}).repsPDataNames{ndInd,durInd}{firstTrialInd};
                    case 'legTrack'
                        thisTrialOptoTime = legTrackOpto.(paramNames{j}).repsOptoTimes{ndInd,durInd,1}(firstTrialInd);
                        thisTrialPDatName = legTrackOpto.(paramNames{j}).repsPDataNames{ndInd,durInd,1}{firstTrialInd};
                end
            % for other parameters, find trial index that matches
            else
                switch paramType{j}
                    case 'fictrac'
                        thisParamOptoTimes = fictracOpto.(paramNames{j}).repsOptoTimes{ndInd,durInd};
                        thisParamPDatNames = fictracOpto.(paramNames{j}).repsPDataNames{ndInd,durInd};
                    case 'legTrack'
                        thisParamOptoTimes = legTrackOpto.(paramNames{j}).repsOptoTimes{ndInd,durInd,1};
                        thisParamPDatNames = legTrackOpto.(paramNames{j}).repsPDataNames{ndInd,durInd,1};
                end

                % get index of this trial, for this parameter, by matching
                %  opto times and pData name
                thisTrial = find((thisParamOptoTimes == thisTrialOptoTime)'...
                    .* (strcmp(thisParamPDatNames,thisTrialPDatName)));
            end

            % get values to plot for this parameter
            % initialize
            thisVal = [];
            switch paramType{j}
                case 'fictrac'
                    thisValNoSmo = fictracOpto.(paramNames{j}).reps{ndInd,durInd}(thisTrial,:);
                    thisTime = fictracOpto.durTs{durInd};
                    % sample rate
                    sampRate = round(1/median(diff(thisTime)));
                    % moving average smoothing; only if no NaNs (doesn't
                    % work if NaNs)
                    if (any(isnan(thisValNoSmo)))
                        thisVal = thisValNoSmo;
                    else
                        thisVal = moveAvgFilt(thisValNoSmo,sampRate,avgWin);  
                    end
                case 'legTrack'
                    % for leg parameters, convert to matrix, each leg as
                    %  column
                    for k = 1:length(legInd)
                        thisLegVal = legTrackOpto.(paramNames{j}).reps{ndInd,durInd,legInd(k)}(thisTrial,:);
                        thisVal(:,k) = thisLegVal';
                    end
                    thisTime = legTrackOpto.durTs{durInd};
            end
            
            % set this subplot
            axHnd(j) = subplot(numParams, numTrials, (j-1) * numTrials + i);
            
            % plot these values
            plot(thisTime, thisVal);

            % lines for when stimulation was on
            line([0 0],yScale(j,:),'Color','k');
            line([whichDur whichDur],yScale(j,:),'Color','k');

            % if this is a FicTrac parameter, add x axis line
            xScale = xlim;
            line(xScale,[0 0], 'Color', 'k');

            % yScale for this parameter
            ylim(yScale(j,:));

            % if this is first parameter for trial, title the plot with the
            %  trial number
            if (j==1)
                ttlStr = sprintf('Trial # %d', thisTrial);
                title(ttlStr);
            end

            % if this is the first trial, add which parameter name to
            % y-axis label
            if (i == 1)
                ylabel(paramNames{j});
            end

            % if this is the last parameter for the trial, label the x-axis
            % with time
            if (j == numParams)
                xlabel('Time (s)');
            end
        end
        % link the time axes together for this trial
        linkaxes(axHnd, 'x');
    end

    % title over whole figure
    % file name: date
    dateName = dataFileName(1:6);
    flyName = dataFileName(8:12);

    plotTitleStr = sprintf('%s %s, ND = %.1f, Stim dur = %.1f s', ...
        dateName, flyName, whichND, whichDur);

    sgtitle(plotTitleStr);

end