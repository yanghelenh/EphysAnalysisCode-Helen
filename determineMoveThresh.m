% determineMoveThresh.m
%
% Function that identifies times when the fly was locomoting vs. not
%  locomoting using a threshold on a smoothed sum of velocities in all 3
%  ball axes. Plots velocities, gives user option to shift threshold,
%  displays threshold. Operates on individual trials'
%  pData and saves back into pData.
%
% INPUTS:
%   pdFile - fullpath to pData file to determine movement threshold on; if
%       [], brings up gui for user to select file
%   thresh - starting movement threshold; value between 0 and 1
%   minBoutLen - minimum length of a bout of stopping or walking, in
%       seconds; those shorter are merged
%   sigmaVel - standard deviation of Gaussian kernel used to smooth 
%       velocity, in seconds
%
% OUTPUTS:
%   none, but updates pData file with new fictrac struct that includes
%       fields:
%     moveParams - parameters for determining movement
%       thresh - final threshold
%       minBoutLen - minimum bout length of movement or stopping, in
%           seconds
%       sigma - standard deviation of Gaussian kernel used to smooth
%           velocity, in seconds
%     moveLog - logical for each time point as moving (1)/not moving (0);
%       does not consider dropInd; values during those times invalid
%     totSpdNormSmo - normalized and smoothed total speed used to compute
%       moveLog; valudes during dropInd not meaningful
%     moveStartInd - indicies into full fictrac for when movement started;
%       values during dropInd not meaningful
%     moveEndInd - indicies into full fictrac for when movement ended;
%       values during dropInd not meaningful
%     moveStartTimes - movement start times, includes start of trial but
%       eliminates any times in dropInd
%     moveEndTimes - movement end times, includes end of trial but
%       elimintates any times in dropInd
%
% CREATED: 9/8/19 - HHY
%
% UPDATED: 
%   9/10/19 - HHY
%   9/27/20 - HHY - updated for ephysData
%

function determineMoveThresh(pdFile, thresh, minBoutLen, sigmaVel)

    % gui to select file if one is not entered as input
    if (isempty(pdFile))
        disp('Select pData.mat file');
        [pdFileName, pdDir] = uigetfile('*.mat');
        pdFile = [pdDir filesep pdFileName];
    end
    
    % load pData file, fictrac strct
    load(pdFile, 'fictracProc', 'fictracParams');
    
    % if movement/not movement has already been dtermined for this pData
    %  file
    if (isfield(fictracProc, 'moveLog'))
        prompt = ['Movement bouts have already been determined for ' ...
            'this pData. Overwrite? (Y/N) '];
        ui = input(prompt,'s');
        if (~strcmpi(ui, 'Y'))
            % stop running this function. don't overwrite 
            disp('Ending determineMoveThresh. Nothing overwritten');
            return;
        end
    end
        
    % interframe interval for fictrac
    ifi = median(diff(fictracProc.t));
    sampRate = 1/ifi; % sample rate for fictrac
    
    % smooth individual axes speeds, more aggressively than in pData
    %   computation
    % sigma for Gaussian kernel smoothing, in samples
    sigmaSamp = round(sigmaVel * sampRate);
    padLen = 3 * sigmaSamp; % pad length, should be longer than sigma
    % convert to integers
    sigmaSamp = int32(sigmaSamp);
    padLen = int32(padLen);
    
    % forward speed
    fwdSpdSmoPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        abs(fictracProc.fwdVel), padLen, sigmaSamp);
    % convert from python to matlab data format
    fwdSpdSmo = cell2mat(cell(fwdSpdSmoPy.tolist()));
    % convert forward velocity into deg/sec (from mm/sec)
    fwdSpdSmo = fwdSpdSmo .* fictracParams.degPerMM;
    
    % slide speed
    slideSpdSmoPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        abs(fictracProc.slideVel), padLen, sigmaSamp);
    % convert from python to matlab data format
    slideSpdSmo = cell2mat(cell(slideSpdSmoPy.tolist()));
    % convert forward velocity into deg/sec (from mm/sec)
    slideSpdSmo = slideSpdSmo .* fictracParams.degPerMM;    
    
    % yaw speed
    yawSpdSmoPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        fictracProc.yawAngSpd, padLen, sigmaSamp);
    % convert from python to matlab data format
    yawSpdSmo = cell2mat(cell(yawSpdSmoPy.tolist()));
    
    % sum movement in all three axes and normalize to max
    totSpdDeg = fwdSpdSmo + slideSpdSmo + yawSpdSmo;
    
    % remove times when fictrac dropped
    validInd = 1:length(totSpdDeg);
    validInd(fictracProc.dropInd) = [];
    validInd = validInd(sigmaSamp:(end-sigmaSamp));
    t = fictracProc.t;
%     t(fictrac.dropInd) = [];
    
    % max speed in deg/sec; remove sigmaSamp from edges b/c of edge effects
    maxSpd = max(totSpdDeg(validInd));
    
    totSpdNorm = totSpdDeg ./ maxSpd;
    
    % smooth normalized total speed again
    totSpdNormSmoPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        totSpdNorm, padLen, sigmaSamp);
    % convert from python to matlab data format
    totSpdNormSmo = cell2mat(cell(totSpdNormSmoPy.tolist())); 
    
    % get total speed, non-smoothed
    totSpd = abs(fictracProc.fwdVel).* fictracParams.degPerMM + ...
        abs(fictracProc.slideVel) .* fictracParams.degPerMM + ...
        fictracProc.yawAngSpd;
%     totSpd(fictrac.dropInd) = [];
    totSpdNotSmo = totSpd ./ max(totSpd(validInd));
    
    % gui for user to select threshold
    tRange = 50;
    tLims = [0 tRange];
    
    f = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
    
    f.Name = 'Set threshold for movement. Close figure when done';
    
    % positions of subplots
    histPos = [0.05 0.3 0.2 0.6];
    smoTotSpdPos = [0.3 0.55 0.68 0.35];
    totSpdPos = [0.3 0.05 0.68 0.35];
    
    % plot histogram of total speeds
    subplot('Position', histPos);
    hist(totSpdNormSmo(validInd),200);
    % plot threshold as vertical black line
    yl = ylim;
    histLine = line([thresh thresh], [0 yl(2)], 'Color', 'black');
    
    % plot total speed, smoothed
    subplotHandles{1} = subplot('Position', smoTotSpdPos);
    plot(t, totSpdNormSmo, 'b');
    % plot threshold as horizontal black line
    totSmoSpdLine = line([0 t(end)], [thresh,thresh],'Color', 'black');
    ylim([0,1]);
    xlim(tLims);
    xlabel('Time (s)');
    title('Total Speed Smoothed');
    
    % plot total speed, not smoothed but normalized
    subplotHandles{2} = subplot('Position', totSpdPos);
    plot(t, totSpdNotSmo, 'r');
    % plot threshold as horizontal black line
    totSpdLine=line([0 t(end)], [thresh,thresh],'Color', 'black');
    ylim([0,1]);
    xlim(tLims);
    xlabel('Time (s)');
    title('Total Speed, not smoothed');
    
    % slider for scrolling x-axis
    tSlider = uicontrol(f, 'Style', 'slider', 'Position', [20 10 350 20]);
    tSlider.Value = 0;
    tSlider.Callback = {@updateTLim, subplotHandles};
    
    % push button to activate video display
    threshBut = uicontrol(f, 'Style', 'pushbutton', 'Position', ...
        [400 10 50 30]);
    threshBut.Callback = @getThresh;
    
    % wait until figure is closed
    uiwait(f);
    
    % use threshold to get moving/not moving bouts
    moveLogical = totSpdNormSmo > thresh;
    
    % indicies where fly transitions between moving and not moving
    transInd = find(diff(moveLogical)) + 1;
    
    % add edges to transitions, so first and last bouts are included
    if transInd(end) == length(moveLogical)
        transInd = [1 transInd];
    else
        transInd = [1 transInd (length(moveLogical)+1)];
    end
    
    % duration of each bout of movement and not movement
    boutDur = diff(transInd);
    
    % cutoff for min bout length, in samples
    minBoutLenSamp = round(minBoutLen * sampRate);
    
    % deal with bouts that are too short (less than minBoutLen) - merge
    % with other, adjacent short bouts or merge into longer sequence
    whichBout = 1;
    while whichBout<=length(boutDur)
        % bout is too short
        if (boutDur(whichBout) < minBoutLenSamp)
            % index into transInd for start of bout
            boutStartInd = whichBout;
            
            % continue until bouts stop being too short
            for k = (whichBout+1):length(boutDur)
                if (boutDur(k) < minBoutLenSamp)
                    whichBout = k;
                else
                    break;
                end
            end
            
            % index into moveLogical for bout transitions
            boutStartMLInd = transInd(boutStartInd);
            % index into tranInd for end of bout
            boutEndInd = whichBout +1;
            % equivalent for index into moveLogical
            boutEndMLInd = transInd(boutEndInd) - 1;
            
            % is this a moving or not moving bout
            boutMoveLogical = moveLogical(boutStartMLInd);
            
            % multiple short bouts, to be merged
            if (whichBout ~= boutStartInd)
                % assign all of these short bouts to 1 longer bout, of the
                % same type as the first bout
                moveLogical(boutStartMLInd:boutEndMLInd) = boutMoveLogical;
            % one short bout, type change   
            else
                moveLogical(boutStartMLInd:boutEndMLInd) = ~boutMoveLogical;
            end
        end
        whichBout = whichBout + 1;
    end
    
    % plot total speed with moving and not moving portions highlighted,
    % with scrolling
    
    % generate patches
    moveStartInd = find(diff(moveLogical) > 0.9) + 1;
    moveStartTimes = t(moveStartInd);
%     moveStartTimes = t((diff(moveLogical) > 0.9) + 1);
    moveEndInd = find(diff(moveLogical) < -0.9);
    moveEndTimes = t(moveEndInd);

    % add start of trial if fly is moving at start
    if (moveLogical(1))
        moveStartTimes = [t(1) moveStartTimes];
        moveStartInd = [1 moveStartInd];
    end
    % add end of trial if fly is moving at end
    if (moveLogical(end))
        moveEndTimes = [moveEndTimes t(end)];
        moveEndInd = [moveEndInd length(moveEndInd)+1];
    end    
    
    % x coordinates of patches
    patchX = [moveStartTimes; moveEndTimes; moveEndTimes; moveStartTimes];
    % y coordinates of patches
    patchY = repmat([0; 0; 1; 1], 1, length(moveStartTimes));
    
    g = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
    subplotHand{1} = subplot(2,1,1);
    plot(t, totSpdNormSmo, 'b');
    hold on;
    patch(patchX, patchY, 'green', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    ylim([0,1]);
    xlim(tLims);
    xlabel('Time (s)');
    
    subplotHand{2} = subplot(2,1,2);
    plot(t, totSpdNotSmo, 'r');
    hold on;
    patch(patchX, patchY, 'green', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    ylim([0,1]);
    xlim(tLims);
    xlabel('Time (s)');
    
    % slider for scrolling x-axis
    tSlider = uicontrol(g, 'Style', 'slider', 'Position', [20 10 350 20]);
    tSlider.Value = 0;
    tSlider.Callback = {@updateTLim, subplotHand};
    
    % remove any move start or end times that fall within fictrac.dropInd
    deleteInd = [];
    for i = 1:length(moveStartInd)
        if ~isempty(find(fictracProc.dropInd==moveStartInd(i),1))
            deleteInd = [deleteInd i];
        end
    end
    moveStartTimes(deleteInd) = [];
    
    deleteInd = [];
    for i = 1:length(moveEndInd)
        if ~isempty(find(fictracProc.dropInd==moveEndInd(i),1))
            deleteInd = [deleteInd i];
        end
    end
    moveEndTimes(deleteInd) = [];
    
    % save threshold, smoothed total speed, moveLogical into fictrac struct
    fictracParams.moveParams.thresh = thresh;
    fictracParams.moveParams.sigma = sigmaVel;
    fictracParams.moveParams.minBoutLen = minBoutLen;
    fictracProc.moveLog = moveLogical;
    fictracProc.totSpdNormSmo = totSpdNormSmo;
    fictracProc.moveStartInd = moveStartInd;
    fictracProc.moveEndInd = moveEndInd;
    
    % note that moveStartTimes and moveEndTimes account for dropInd and
    %  ignore any times that fall during this time period
    fictracProc.moveStartTimes = moveStartTimes;
    fictracProc.moveEndTimes = moveEndTimes;
    
    % update pData
    save(pdFile, 'fictrac', 'fictracParams', '-append');
    
    
% HELPER FUNCTIONS FOR PLOTTING
    % function for scroll bar in threshold selection figure
    function updateTLim(src, event, subplotHandles)
        for j = 1:length(subplotHandles)
            xlim(subplotHandles{j}, ...
                [tSlider.Value * (t(end)-tRange),...
                tSlider.Value * (t(end)-tRange) + tRange]);
            ylim(subplotHandles{j}, [0,1]);
        end
    end

    % function to get user to pick a threshold
    function getThresh(src, event)
        DONE = 0;
        
        while ~DONE
            [xVal, ~] = ginput(1);

            if ~isempty(xVal)
                thresh = xVal;
            else
                DONE = 1;
            end
            histLine.XData = [thresh, thresh];
            totSmoSpdLine.YData = [thresh, thresh];
            totSpdLine.YData = [thresh,thresh];
        end            
    end

end