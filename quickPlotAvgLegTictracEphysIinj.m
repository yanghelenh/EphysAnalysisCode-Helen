% script for quick plots of avgLegFictracEphysIinj

figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 3;

lgdStr = {};

for i = 1:length(amps)
    % moving average smooth FicTrac data
    movAvgFwdVel = moveAvgFilt(fictracIinj.fwdVel.means{i,durInd},sampRate,avgWin);
    plot(fictracIinj.durTs{durInd},movAvgFwdVel);
    lgdStr{i} = sprintf('Amp = %d pA',amps(i));
end

yScale = [-5 20];
xScale = [fictracIinj.durTs{durInd}(1) fictracIinj.durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('FwdVel (mm/s)');

ttlStr = sprintf('Forward Velocity, %.1f sec stim', durs(durInd));

title(ttlStr);

legend(lgdStr);

%%
figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 3;

lgdStr = {};

for i = 1:length(amps)
    % moving average smooth FicTrac data
    movAvgYawVel = moveAvgFilt(fictracIinj.yawAngVel.means{i,durInd},sampRate,avgWin);
    plot(fictracIinj.durTs{durInd},movAvgYawVel);
    lgdStr{i} = sprintf('Amp = %d pA',amps(i));
end

yScale = [-200 200];
xScale = [fictracIinj.durTs{durInd}(1) fictracIinj.durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('YawAngVel (deg/s)');

ttlStr = sprintf('Yaw Velocity, %.1f sec stim', durs(durInd));

title(ttlStr);

legend(lgdStr);

%%
figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 2;

lgdStr = {};

for i = 1:length(amps)
    % moving average smooth FicTrac data
    movAvgYawVel = moveAvgFilt(fictracIinj.slideVel.means{i,durInd},sampRate,avgWin);
    plot(fictracIinj.durTs{durInd},movAvgYawVel);
    lgdStr{i} = sprintf('Amp = %d pA',amps(i));
end

yScale = [-5 5];
xScale = [fictracIinj.durTs{durInd}(1) fictracIinj.durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('SlideVel (mm/s)');

ttlStr = sprintf('Side Velocity, %.1f sec stim', durs(durInd));

title(ttlStr);

legend(lgdStr);