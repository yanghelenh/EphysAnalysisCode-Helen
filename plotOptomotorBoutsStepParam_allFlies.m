% plotOptomotorBoutsStepParam_allFlies.m
%
% Function that plots the output of
%  extractOptomotorBoutsLegStepParams_fly() for one, user specified step
%  parameter. User can select subset of bouts based on peak velocity or
%  change in velocity.
%
% INPUTS:
%   condPeakVel - struct of parameters for selecting bouts based on peak
%       velocity (smo), [] for no conditions
%     whichParam - cell array of combo of 'yaw', 'fwd', 'slide'
%     cond - condition for each parameter, matched
%   condChangeVel - struct of parameters for selecting bouts based on the
%       change in velocity (smo) between the peak and start, [] for no
%       conditions
%     whichParam - cell array of combo of 'yaw', 'fwd', 'slide'
%     cond - condition for each parameter, matched
%   datDir - directory with output files
%   whichParam - which step parameter to plot
%   whichPhase - which step phase to plot
%   whichVel - which optomotor velocity condition to plot
%   NDs - which NDs to plot
%   minNumBouts - minimum number of bouts for the fly to be included
%   yScale - scale for plots, as [min max]
%   plotIndiv - boolean for whether to plot individual flies
%   plotDiff - boolean for whether to plot values as difference from 1st
%       step
%   plotRelNoStim - boolean for whether to plot values as difference from
%       no optogenetic stimulation (ND = -1)
%   numSteps - number of steps relative to yaw peak to plot (on either
%       side)
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 4/17/24 - HHY
%
% UPDATED:
%   4/17/24 - HHY
%
function plotOptomotorBoutsStepParam_allFlies(condPeakVel, ...
    condChangeVel, datDir, whichParam, whichPhase, whichVel, NDs, ...
    minNumBouts, yScale, plotIndiv, plotDiff, plotRelNoStim, numSteps)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % circular step parameters
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveLegStepParamByCond_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select Step Param files', datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    numNDs = length(NDs);

    % initialize
    allFliesMeans = zeros(numSteps * 2+1, 6, numNDs, numFlies);
    allFliesSEMs = zeros(numSteps* 2+1, 6, numNDs, numFlies);
    allFliesStd = zeros(numSteps* 2+1, 6, numNDs, numFlies);
    allFliesN = zeros(numSteps* 2+1, 6, numNDs, numFlies);

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        if (strcmpi(whichPhase,'swing'))
            load(outputFullPath, 'selSwingParams', 'selOptoVels', ...
                'selNDs', 'boutPeakVelSmo', 'boutStartVelSmo', ...
                'numBouts', 'maxNumSteps');
            thisParamAll = selSwingParams.(whichParam);
        else
            load(outputFullPath, 'selStanceParams', 'selOptoVels', ...
                'selNDs', 'boutPeakVelSmo', 'boutStartVelSmo', ...
                'numBouts','maxNumSteps');
            thisParamAll = selStanceParams.(whichParam);
        end 

        % conditioning on peak velocity
        if (~isempty(condPeakVel))
            peakCondLog = true(numBouts,1);
            for j = 1:length(condPeakVel.whichParam)
                thisCondParam = boutPeakVelSmo.(condPeakVel.whichParam{j});

                thisCondLog = eval(['thisCondParam' condPeakVel.cond{j}]);

                % combine with other conditions
                peakCondLog = peakCondLog & thisCondLog;
            end
        else
            peakCondLog = true(numBouts,1);
        end

        % conditioning on change in velocity
        if (~isempty(condChangeVel))
            changeCondLog = true(numBouts,1);
            for j = 1:length(condChangeVel.whichParam)
                thisCondParam = boutPeakVelSmo.(condChangeVel.whichParam{j}) - ...
                    boutStartVelSmo.(condChangeVel.whichParam{j});

                thisCondLog = eval(['thisCondParam' condChangeVel.cond{j}]);

                % combine with other conditions
                changeCondLog = changeCondLog & thisCondLog;
            end
        else
            changeCondLog = true(numBouts,1);
        end

        % select valid optomotor velocity
        optoVelLog = selOptoVels == whichVel;

        % get logical for all
        allLog = peakCondLog & changeCondLog & optoVelLog;

        % apply logical to param, NDs
        thisParamSel = thisParamAll(:,:,allLog);
        thisNDsSel = selNDs(allLog);

        % get indices of steps to consider
        pkStepInd = maxNumSteps + 1;
        stepInds = (pkStepInd-numSteps):(pkStepInd+numSteps);

        % if plotting relative to no stim, get values for that condition
        if plotRelNoStim
            noStimMeans = zeros(numSteps * 2+1, 6);
            for j = 1:length(stepInds)
                thisStepInd = stepInds(j);
    
                for k = 1:6
                    thisNDLog = (thisNDsSel == -1);
    
                    thisTandLegs = thisParamSel(thisStepInd, k, thisNDLog);

                    % get difference if needed
                    if plotDiff
                        thisTandLegs = thisTandLegs - thisParamSel(stepInds(1),k,thisNDLog);
                    end

                    % remove NaNs
                    thisTandLegs(isnan(thisTandLegs)) = [];

                    thisN = length(rmoutliers(squeeze(thisTandLegs)));

                    if (thisN >= minNumBouts)
    
                        % outlier removal before calculating mean and std
                        % get mean and std dev, account for if circular 
                        if(any(strcmpi(whichParam, circStepParams)))
                            % convert this parameter to radians
                            thisTandLegs = deg2rad(thisTandLegs);
                            % get circular mean
                            thisMean = rad2deg(circ_mean(rmoutliers(squeeze(thisTandLegs))));
                        % if not circular, compute regular mean and std
                        else
                            thisMean = mean(rmoutliers(squeeze(thisTandLegs)));
                        end

                    else
                        thisMean = nan;
                    end

                    % save into output
                    noStimMeans(j,k) = thisMean;
                end
            end
        end


        % loop through all NDs, get values for this fly, this ND
        for j = 1:numNDs
            thisNDLog = thisNDsSel == NDs(j);

            thisParamND = thisParamSel(:,:,thisNDLog);

            % loop through all time points
            for k = 1:length(stepInds)
                thisStepInd = stepInds(k);
                
                % loop through all legs
                for l=1:6
                    thisTandLegs = thisParamND(thisStepInd,l,:);

                    % get difference if needed
                    if plotDiff
                        thisTandLegs = thisTandLegs - thisParamND(stepInds(1),l,:);
                    end
    
                    % remove NaNs
                    thisTandLegs(isnan(thisTandLegs)) = [];

                    % n for this time point and legs
                    thisN = length(rmoutliers(squeeze(thisTandLegs)));

                    % check that there are enough data points for inclusion
                    if (thisN>=minNumBouts)
    
                        % outlier removal before calculating mean and std
                        % get mean and std dev, account for if circular 
                        if(any(strcmpi(whichParam, circStepParams)))
                            % convert this parameter to radians
                            thisTandLegs = deg2rad(thisTandLegs);
                            % get circular mean
                            thisMean = rad2deg(circ_mean(rmoutliers(squeeze(thisTandLegs))));
                            % get circular std
                            thisStd = rad2deg(circ_std(rmoutliers(squeeze(thisTandLegs))));
                        % if not circular, compute regular mean and std
                        else
                            thisMean = mean(rmoutliers(squeeze(thisTandLegs)));
                            thisStd = std(rmoutliers(squeeze(thisTandLegs)));
                        end
    
                        % get SEM
                        thisSEM = thisStd / sqrt(thisN);
    
                        % add to output matrix
                        if plotRelNoStim
                            allFliesMeans(k,l,j,i) = thisMean - noStimMeans(k,l);
                        else
                            allFliesMeans(k,l,j,i) = thisMean;
                        end
                        allFliesSEMs(k,l,j,i) = thisSEM;
                        allFliesStd(k,l,j,i) = thisStd;
                        allFliesN(k,l,j,i) = thisN;

                    % if there aren't enough data points
                    else
                        allFliesMeans(k,l,j,i) = nan;
                        allFliesSEMs(k,l,j,i) = nan;
                        allFliesStd(k,l,j,i) = nan;
                        allFliesN(k,l,j,i) = thisN;
                    end
                end
            end
        end
    end

    % compute mean and SEM across flies
    totMeans = zeros(numSteps * 2+1, 6, numNDs);
    totSEMs = zeros(numSteps* 2+1, 6, numNDs);
    totFliesN = zeros(numSteps* 2+1, 6, numNDs);
    % loop through all steps
    for i = 1:length(stepInds)
        % loop through all legs
        for j = 1:6
            % loop through all NDs
            for k = 1:numNDs
                % values for all flies for this step, leg, ND
                thisStepLegND = allFliesMeans(i,j,k,:);

                thisStepLegND(isnan(thisStepLegND)) = [];

                thisStepLegND = squeeze(thisStepLegND);

                % compute mean and std
                % if circular
                if(any(strcmpi(whichParam, circStepParams)))
                    thisStepLegND = deg2rad(thisStepLegND);

                    % compute circular mean and std
                    thisMean = rad2deg(circ_mean(thisStepLegND));
                    thisStd = rad2deg(circ_std(thisStepLegND));

                % not circular    
                else
                    thisMean = mean(thisStepLegND);
                    thisStd = std(thisStepLegND);
                end

                % compute N and SEM
                thisN = length(thisStepLegND);
                thisSEM = thisStd / sqrt(thisN);

                % add to output matricies
                totMeans(i,j,k) = thisMean;
                totSEMs(i,j,k) = thisSEM;
                totFliesN(i,j,k) = thisN;
            end
        end
    end

    % get step time points (for x-axis): 0 for step at peak
    stepTPts = [fliplr((1:numSteps) * -1), 0, 1:numSteps];
    

    % initialize figure
    f = figure;
    c = colormap('lines');

    % legend, specifies ND
    legendStr = cell(numNDs,1);
    for i = 1:numNDs
        legendStr{i} = sprintf('ND=%.1f',NDs(i));
    end

    % handles to plot of mean across flies
    meanHdl = zeros(numNDs,1);

    for i = 1:6
        subplot(3,2,subInd(i));
        hold on;

        % loop through all NDs
        for j = 1:numNDs
            % plot mean across flies
            meanHdl(j) = errorbar(stepTPts', totMeans(:,i,j), ...
                totSEMs(:,i,j), 'Marker', 'x','LineWidth',2, ...
                'CapSize', 0, 'Color', c(j,:));

            % plot mean for each individual fly (no error bar)
            if plotIndiv
                plot(stepTPts', squeeze(allFliesMeans(:,i,j,:)),'LineWidth',1,...
                    'Color', c(j,:));
            end

            hold on;
        end


        % axis scale and label

        % for AEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-0.9 -0.5]);
%         elseif (subInd(i)==3 || subInd(i)==4)
%             ylim([-0.3 0.1]);
%         else
%             ylim([0.15 0.55]);
%         end

        % for PEP
%         if (subInd(i)==1 || subInd(i)==2)
%             ylim([-0.6 -0.2]);
%         elseif (subInd(i)==3 || subInd(i)==4)
%             ylim([0 0.4]);
%         else
%             ylim([0.45 0.85]);
%         end

        ylim(yScale);

        xScale = xlim;
        xScale(1) = xScale(1) - (0.1 * (stepTPts(end)-stepTPts(1)));
        xScale(2) = xScale(2) + (0.1 * (stepTPts(end)-stepTPts(1)));
        xlim(xScale);

        xticks(stepTPts);

        % add line at y = 0 if plotting difference
        if (plotDiff)
            line(xScale,[0,0],'Color','k', 'LineWidth', 1);
        end


        xlabel('Steps');
        ylabel(whichParam);

        % legend
        legend(meanHdl,legendStr);
    end

    sgtitle(sprintf('%s, %s, vel=%d', whichParam, whichPhase, whichVel));
end