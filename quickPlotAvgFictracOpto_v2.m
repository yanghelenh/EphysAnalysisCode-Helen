datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData/220520_LegFictracOpto_1Fly';    
% prompt user to select output files from saveLegStepParamCond_indpt()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select avgFicTrac files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

NDInd = 3;

figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 1;

for i = 1:numFiles
    
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];


        % load data
        load(outputFullPath);

        plot(fictracOpto.durTs{durInd},fictracOpto.fwdVel.means{NDInd,durInd});

        hold on;
end

yScale = [-5 20];
xScale = [fictracOpto.durTs{durInd}(1) fictracOpto.durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('FwdVel (mm/s)');

ttlStr = sprintf('Forward Velocity, %.1f sec stim, %.1f ND', durs(durInd), NDs(NDInd));

title(ttlStr);

%%

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData/220520_LegFictracOpto_1Fly';    
% prompt user to select output files from saveLegStepParamCond_indpt()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select avgFicTrac files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

NDInd = 3;

figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 2;

for i = 1:numFiles
    
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];


        % load data
        load(outputFullPath);

        plot(fictracOpto.durTs{durInd},fictracOpto.yawAngVel.means{NDInd,durInd});

        hold on;
end

yScale = [-200 200];
xScale = [fictracOpto.durTs{durInd}(1) fictracOpto.durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('YawVel (deg/s)');

ttlStr = sprintf('Yaw Velocity, %.1f sec stim, %.1f ND', durs(durInd), NDs(NDInd));

title(ttlStr);