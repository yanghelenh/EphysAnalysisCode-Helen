% script for quick plots of avgFictracOpto

figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 2;

for i = 1:length(NDs)
    % moving average smooth FicTrac data
    movAvgFwdVel = moveAvgFilt(fwdVelData.means{i,durInd},sampRate,avgWin);
    plot(durTs{durInd},movAvgFwdVel);
end

yScale = [-5 20];
xScale = [durTs{durInd}(1) durTs{durInd}(end)];
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

% legend({'ND = 3', 'ND = 2', 'ND = 1.3','ND = 0.6', 'ND = 0'});
legend({'ND = 2.0', 'ND = 1.3','ND = 0.6', 'ND = 0'});
% legend('ND = 0');
% legend({'ND = 3', 'ND = 2'});
% legend({'ND = 1.3','ND = 0.6', 'ND = 0'});
% legend({'ND = 0.6', 'ND = 0'});
% legend('ND = 0.6');
% legend({'ND = 3', 'ND = 2', 'ND = 1.3','ND = 0.6'});

%%

figure;
hold on;

sampRate = 1000;
avgWin = 0.1;

durInd = 2;

for i = 1:length(NDs)
    % moving average smooth FicTrac data
    movAvgYawVel = moveAvgFilt(yawVelData.means{i,durInd},sampRate,avgWin);
    plot(durTs{durInd},movAvgYawVel);
end

yScale = [-200 200];
xScale = [durTs{durInd}(1) durTs{durInd}(end)];
ylim(yScale);
xlim(xScale);

% stimulus lines
line([0 0],yScale,'Color','k');
line([durs(durInd) durs(durInd)],yScale,'Color','k');

% x-axis line
line(xScale,[0 0], 'Color', 'k');

xlabel('Time (s)');
ylabel('YawVel (deg/s)');

ttlStr = sprintf('Yaw Velocity, %.1f sec stim', durs(durInd));

title(ttlStr);

% legend({'ND = 3', 'ND = 2', 'ND = 1.3','ND = 0.6', 'ND = 0'});
legend({'ND = 2.0', 'ND = 1.3','ND = 0.6', 'ND = 0'});
% legend({'ND = 3', 'ND = 2'});
% legend('ND = 0');
% legend({'ND = 1.3','ND = 0.6', 'ND = 0'});
% legend({'ND = 0.6', 'ND = 0'});
% legend('ND = 0.6');
% legend({'ND = 3', 'ND = 2', 'ND = 1.3','ND = 0.6'});
