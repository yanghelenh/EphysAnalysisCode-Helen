% interactEphysGetSpikes.m
%
% Function that detects spikes, computes spike rate, and computes
%  median-filtered voltage trace on each preprocessed ephys trace. Uses GUI
%  for user to select parameters for spike calling. Runs on selected pData 
%  files. Must have run preprocess() before running this.
%
% Replaces ephysGetSpikes_all.m
%
% INPUTS:
%   none, but prompts user to select pData files and spike-calling
%       parameters
%
% OUTPUTS:
%   none, but saves spike info back into same pData file
%
% CREATED: 6/29/22 - HHY
%
% UPDATED:
%   6/30/22 - HHY
%
function interactEphysGetSpikes()

    % initial parameter values
    ephysSpikes.params.dvdtThresh = 3000; 
    % 2 ms min time b/w spikes
    ephysSpikes.params.refractoryPeriod = 2/1000; 
    % 50 ms median filter to remove spikes
    ephysSpikes.params.medFiltOrder = 0.05; 

    % ranges for parameters
    ephysSpikes.paramRanges.dvdtThresh = [0 15000];
    ephysSpikes.paramRanges.refractoryPeriod = [0 10/1000]; % in sec
    ephysSpikes.paramRanges.medFiltOrder = [0 0.5]; % in sec

    % save initial values
    ephysSpikesInit = ephysSpikes;

    paramNames = fieldnames(ephysSpikes.params);

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
    
    for i = 1:numPDataFiles

        % initialize logical for done button press
        doneButtonPushed = 0;
        
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];
        
        % load pData
        load(pDataFullPath, 'exptCond');
        
        % check that pData has ephys data, otherwise, skip
        if (contains(exptCond, 'Ephys', 'IgnoreCase', true))
            load(pDataFullPath, 'ephysData');

            % reset ephysSpikes parameter values
            ephysSpikes = ephysSpikesInit;

            % some parameters for interactive figure
            tRange = 1; % in sec, amount of data to display
            xMax = ephysData.t(end); % max value for t

            % slider parameters
            sldXPos = 1200;
            sldYPosStart = 850;
            sldYPosEnd = 50;
            sldHeight = 20;
            sldWidth = 300;
            numSld = length(fieldnames(ephysSpikes.params));
            sldYSpace = round((sldYPosStart - sldYPosEnd) / numSld);
            sldTMajorStep = (tRange/xMax) * 0.9;
            sldTMinorStep = 0.1 * sldTMajorStep;

            % get ylimits - on scaled voltage; 5% buffer
            sVMin = min(ephysData.scaledVoltage);
            sVMax = max(ephysData.scaledVoltage);
            yMinLim = sVMin - (0.05 * abs(sVMax - sVMin));
            yMaxLim = sVMax + (0.05 * abs(sVMax - sVMin));

            % GET SPIKES %

            % get spike calls, using initial parameter values
            spikeStartInds = detectSpikes(...
                ephysData.scaledVoltage, ephysData.t, ...
                ephysSpikes.params.dvdtThresh, ...
                ephysSpikes.params.refractoryPeriod);
            % turn spike start inds into spike times and voltage values
            spikeStartTs = ephysData.t(spikeStartInds);
            spikeStartVm = ephysData.scaledVoltage(spikeStartInds);

            % get median-filtered voltage trace, using initial parameter
            %  values
            medFiltV = filtEphysNoSpikes(...
                ephysData.scaledVoltage, ephysData.t, ...
                ephysSpikes.params.medFiltOrder);


             % initialize figure
            f = figure('Position', [20 20 1600 920]);
            
            % plot scaled voltage with spike starts as red x on top of
            scaledVoltAx = subplot('Position', [0.05 0.6 0.6 0.3]);
            % scaled voltage
            plot(ephysData.t, ephysData.scaledVoltage);
            hold on;
            % spike starts
            plot(spikeStartTs, spikeStartVm, 'x', 'LineStyle','none');

            % scale and label
            xlim([0 tRange]);
            ylim([yMinLim yMaxLim]);
            title('Scaled Voltage')
            xlabel('Time (s)');
            ylabel('mV');
            
            % plot median-filed voltage trace
            medFiltVAx = subplot('Position', [0.05 0.15 0.6 0.3]);
            plot(ephysData.t, medFiltV);
            xlim([0 tRange]);
            ylim([yMinLim yMaxLim]);
            title('Median Filtered Voltage');
            xlabel('Time (s)');
            ylabel('mV');

            % lock x axes of 2 plots together
            linkaxes([scaledVoltAx medFiltVAx],'x');

            
            % Set up sliders
            % text for names of parameter sliders
            txtAx = axes('Position', [0 0 1 1], 'Visible', 'off');
            set(gcf, 'CurrentAxes', txtAx);
            
            % text for slider labels
            for j = 1:length(paramNames)
                thisTxtPos = sldYPosStart - sldYSpace * (j-1) + 10;
                
                text(sldXPos - 120, thisTxtPos, paramNames{j}, ...
                    'Units', 'pixels', 'FontSize', 12);
            end
            
            % text for values of parameter sliders
            allTxtH = {}; % handles to text obj
            for j = 1:length(paramNames)
                thisTxtPos = sldYPosStart - sldYSpace * (j-1) + 10;
                
                thisDispVal = num2str(ephysSpikes.params.(paramNames{j}));
                
                allTxtH{j} = text(sldXPos + 320, thisTxtPos, thisDispVal, ...
                    'Units', 'pixels', 'FontSize', 12);
            end
            
            % slider adjusting x-axis view
            tSlider = uicontrol(f, 'Style', 'slider', 'Position', ...
                [100 50 600 20]);
            tSlider.Value = 0;
            tSlider.SliderStep = [sldTMinorStep sldTMajorStep];
            tSlider.Callback = @updateTLim;
             
    
                
            % sliders for adjusting each of the parameter values
            % initialize cell array for slider objects
            allSld = {};
            % loop through all parameters, creating sliders
            for j = 1:length(paramNames)
                allSld{j} = uicontrol(f, 'Style', 'slider');
    
                thisSldYPos = sldYPosStart - sldYSpace * (j - 1);
    
                allSld{j}.Position = [sldXPos thisSldYPos sldWidth sldHeight];
    
                allSld{j}.Value = ephysSpikes.params.(paramNames{j});
                thisParamRange = ephysSpikes.paramRanges.(paramNames{j});
                allSld{j}.Max = thisParamRange(2);
                allSld{j}.Min = thisParamRange(1);

                allSld{j}.Callback = {@updateGraph, j};
            end

            % button for when done
            doneButton= uicontrol(f, 'Style', 'pushbutton', 'String', ...
                'Done', 'Position', [100 30 50 20]);
            doneButton.Callback = @donePushed;

            % loop until user hits next leg/done button
            while ~doneButtonPushed
                pause(0.1);
            end

            % compute spike rate
            ephysSpikes.spikeRate = computeSpikeRate(...
                spikeStartInds, ephysData.t);
           
            % save into ephysSpikes struct
            ephysSpikes.startInd = spikeStartInds;
            ephysSpikes.medFiltV = medFiltV;

            % copy over time vector
            ephysSpikes.t = ephysData.t;

            % close figure
            close(f);
            
            % save back into pData
            save(pDataFullPath, 'ephysSpikes', '-append');
            
            % display
            fprintf('Saved ephysSpikes for %s!\n', pDataName);
            
        else
            % display
            fprintf('%s does not have ephys data\n', pDataName);
        end
    end

        % functions for figure, down here b/c can't be in loop
    
    % function to update display region for plot, every time that
    %  slider is moved
    function updateTLim(src, event)
        xlim(medFiltVAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
        ylim(medFiltVAx, [yMinLim yMaxLim]);
    end

    % function for updating spike calls, figure and parameters
    function updateGraph(src, event, nameInd)

        % get parameter value
        thisVal = allSld{nameInd}.Value;

        % change the apppropriate parameter value in thisLegRevParams
        ephysSpikes.params.(paramNames{nameInd}) = thisVal;

        % get spike calls, using new parameter values
        spikeStartInds = detectSpikes(...
            ephysData.scaledVoltage, ephysData.t, ...
            ephysSpikes.params.dvdtThresh, ...
            ephysSpikes.params.refractoryPeriod);
        % turn spike start inds into spike times and voltage values
        spikeStartTs = ephysData.t(spikeStartInds);
        spikeStartVm = ephysData.scaledVoltage(spikeStartInds);

        % get new median filtered V
        medFiltV = filtEphysNoSpikes(...
            ephysData.scaledVoltage, ephysData.t, ...
            ephysSpikes.params.medFiltOrder);

%         % get current xlims
%         currXLims = xlim(scaledVoltAx);

        % clear scaled voltage axes
        cla(scaledVoltAx);

        % plot
        plot(scaledVoltAx,ephysData.t,ephysData.scaledVoltage);
        hold on;
        plot(scaledVoltAx,spikeStartTs,spikeStartVm,'x','LineStyle','none');

        ylim(scaledVoltAx, [yMinLim yMaxLim]);

% 
%         % scale and label
%         xlim(currXLims);
%         title('Scaled Voltage')
%         xlabel('Time (s)');
%         ylabel('mV');

        % clear med filt axes
        cla(medFiltVAx);

        % plot
        plot(medFiltVAx,ephysData.t, medFiltV);

%         % scale and label
%         xlim(currXLims);
%         title('Median Filtered Voltage');
%         xlabel('Time (s)');
%         ylabel('mV');

%         % lock x axes of 2 plots together
%         linkaxes([scaledVoltAx medFiltVAx],'x');

        xlim(medFiltVAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
        ylim(medFiltVAx, [yMinLim yMaxLim]);

        % update display value around slider
        allTxtH{nameInd}.String = num2str(thisVal);
    end

    % function when done button pushed
    function donePushed(src, event)
        % toggle logical to done, stops loop
        doneButtonPushed = 1;
    end
end