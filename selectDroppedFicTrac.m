% selectDroppedFicTrac.m
%
% Function to manually select times in experiment where FicTrac was not
%  tracking the ball. Loads FicTrac data from pData file. Saves back into
%  same pData file.
%
% Adapted from function of same name in 2PAnalysisCode-Helen
%
% INPUT: 
%   none, but user selects trial folder to process
%
% OUTPUT:
%   none, but saves dropInd variable into pData file
%
% Created: 9/11/20 - HHY
%
% Updated: 
%   9/14/20 - HHY
%

function selectDroppedFicTrac()

    T_RANGE = 50; % how many seconds of FicTrac data to display at once

    % ask user to select pData file
    [fileName, filePath] = uigetfile('*.mat', 'Select pData file', ...
        pDataDir());
    
    % full path to pData file
    fullPath = [filePath fileName];
    
    % load pData file
    load(fullPath, 'exptCond');
    
    % if this pData file describes experiment with FicTrac data
    if (contains(exptCond, 'Fictrac', 'IgnoreCase', true))
        % load fictrac struct only if experiment has it
        load(fullPath, 'fictrac');
        
        % check whether dropped frames previously selected (i.e. dropInd
        %  field of fictrac struct exists); if so, prompt user to ask
        %  whether to overwrite
        if (isfield(fictrac, 'dropInd'))
            prompt = ['FicTrac dropping times have already been '...
                'selected for this trial.' ' Overwrite? (Y/N) '];
            ui = input(prompt,'s');
            if (~strcmpi(ui, 'Y')) 
                % stop running this function. don't overwrite 
                disp('Ending selectDroppedFicTrac. Nothing overwritten');
                return;
            end
        end
        
        % pre-generate dropInd variable and save into fictrac struct of
        %  pData
        fictrac.dropInd = [];
        save(fullPath, 'fictrac', '-append');
        
        % generate plot for selecting times of FicTrac dropping
        
        % determine starting xLim on all plots
        tLims = [0 T_RANGE];
        
        % end time of fictrac data
        xMax = fictrac.t(end);

        f = figure;
        hold on;
        
        % plot fwd velocity
        fwdVelSubplot = subplot(2, 1, 1);
        fwdVelPlot = plot(fictrac.t, fictrac.fwdVel);
        xlim(tLims);
%         ylim([-10 30]);
        xlabel('Time (sec)');
        ylabel('mm/sec');
        title('Forward Velocity');
        
        % plot yaw angular velocity
        yawVelSubplot = subplot(2, 1, 2);
        yawVelPlot = plot(fictrac.t, fictrac.yawAngVel);
%         ylim([-500 500]);
        xlim(tLims);
        xlabel('Time (sec)');
        ylabel('deg/sec');
        title('Rotational Velocity');
        
        % link x-axes
        linkaxes([fwdVelSubplot, yawVelSubplot], 'x');
        
        % create brush object
        brushobj = brush(f);
        brushobj.Enable = 'off'; % start with brushobj off
        
        % slider for scrolling x-axis
        tSlider = uicontrol(f, 'Style', 'slider', 'Position', ...
            [20 10 400 20]);
        tSlider.Value = 0;
        tSlider.Callback = @updateTLim;

        % push button to activate/inactivate brushing to select dropped
        %  times
        brushBut = uicontrol(f, 'Style', 'togglebutton', 'Position', ...
            [500 10 50 30]);
        brushBut.Callback = @brushButtonPush;
        brushBut.Value = brushBut.Min; % start with button deselected
        
    else
        disp('The selected pData file does not contain FicTrac data');
        return;
        
    end

    % functions for interacting with plot
    % function to change time axis with slider
    function updateTLim(src, event)
        xlim(fwdVelSubplot, [tSlider.Value * (xMax-T_RANGE),...
                tSlider.Value * (xMax-T_RANGE) + T_RANGE]);
        xlim(yawVelSubplot, [tSlider.Value * (xMax-T_RANGE),...
                tSlider.Value * (xMax-T_RANGE) + T_RANGE]); 
    end

    % function to start/stop brushing, save values
    function brushButtonPush(src, event)
        % if toggle button is selected
        if (brushBut.Value == brushBut.Max)
            % enable brushing
            brushobj.Enable = 'on';
        
        % if toggle button is deselected
        elseif (brushBut.Value == brushBut.Min)
            % disable brushing
            brushobj.Enable = 'off';
            
            % get brushed values
            brushFwdVelInd = find(get(fwdVelPlot, 'BrushData'));
            brushYawVelInd = find(get(yawVelPlot, 'BrushData'));
            
            % FicTrac dropped values are union of values selected on 2
            %  plots
            fictrac.dropInd = union(brushFwdVelInd, brushYawVelInd);
            
            % save into pData file
            save(fullPath, 'fictrac', '-append');
        end
    end
end