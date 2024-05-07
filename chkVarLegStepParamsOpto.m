% chkVarLegStepParamsOpto.m
%
% Quick function to check the variance of a given leg step parameter, for
%  one or more flies. Input is output of extractLegStepParamsOpto_fly() or
%  extractOptomotorLegStepParamsOptoCond_fly()
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter
%   whichPhase - which step phase
%   whichVelDur - which velocity or duration condition
%   whichND - which ND
%
% OUTPUTS:
%   allFliesStdDev - matrix of numFlies x 6 legs for std dev of parameter
%       values, for each fly and leg
%
% CREATED: 5/7/24 - HHY
%
% UPDATED:
%   5/7/24 - HHY
%
function allFliesStdDev = chkVarLegStepParamsOpto(datDir, whichParam, ...
    whichPhase, whichVelDur, whichND)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files 
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select Step Param files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % preallocate
    allFliesStdDev = zeros(numFlies,6); % 6 for number of legs

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        flyName = outName(1:19); % name of fly

        % get variable from file
        theseVars = who('-file', outputFullPath);

        if (any(strcmpi(theseVars,'condKeyVels')))
            % load data
            load(outputFullPath, 'legStepsOptoAll', ...
                'legStepsOptoStdDev', 'condKeyVels', 'condKeyNDs');
            condKeyVD = condKeyVels';
        else
            load(outputFullPath, 'legStepsOptoAll', ...
                'legStepsOptoStdDev', 'condKeyDurs', 'condKeyNDs');
            condKeyVD = condKeyDurs;
        end

        % get index of particular condition
        thisCondLog = (condKeyVD == whichVelDur) & ...
            (condKeyNDs == whichND);

        if (strcmpi(whichPhase, 'stance'))
            if (any(strcmpi(theseVars,'condKeyVels')))
                selSteps = (legStepsOptoAll.stance.visVel == whichVelDur) & ...
                    (legStepsOptoAll.stance.optoND == whichND);
            else
                selSteps = (legStepsOptoAll.stance.optoDur == whichVelDur) & ...
                    (legStepsOptoAll.stance.optoND == whichND);
            end
            thisWhichLeg = legStepsOptoAll.stance.stepWhichLeg(selSteps);
            thisParam = legStepsOptoAll.stance.(whichParam)(selSteps);
            thisStdDev = legStepsOptoStdDev.stance.(whichParam)(thisCondLog,:);

        elseif (strcmpi(whichPhase, 'swing'))
            if (any(strcmpi(theseVars,'condKeyVels')))
                selSteps = (legStepsOptoAll.swing.visVel == whichVelDur) & ...
                    (legStepsOptoAll.swing.optoND == whichND);
            else
                selSteps = (legStepsOptoAll.swing.optoDur == whichVelDur) & ...
                    (legStepsOptoAll.swing.optoND == whichND);
            end
            thisWhichLeg = legStepsOptoAll.swing.stepWhichLeg(selSteps);
            thisParam = legStepsOptoAll.swing.(whichParam)(selSteps);
            thisStdDev = legStepsOptoStdDev.swing.(whichParam)(thisCondLog,:);
        end

        % add std dev to tracker across flies
        allFliesStdDev(i,:) = thisStdDev;

        % plot histograms for this fly
        figure;

        % loop across all legs, one subplot per leg
        for j = 1:6
            subplot(3,2,subInd(j));
            hold on;

            thisLegVals = thisParam(thisWhichLeg == j);
            histogram(thisLegVals,20);
        end

        sgtitle(sprintf('%s, %s, %s',flyName, whichParam, whichPhase));
    end
end