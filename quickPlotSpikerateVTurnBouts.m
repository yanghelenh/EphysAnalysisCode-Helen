% quickPlotSpikerateVTurnBouts.m
%
% 7/18/23 - Quick and dirty plotting function relating spike rate to
% turning bouts. Plots output of saveSpikerate_bouts.m

for i = 1:5

    figure;
    
    scatter(allYaw(i,:), allFwd(i,:), [], allSpikerate(i,:), 'filled');
    colormap(gca, 'parula');
    caxis([0 100])
    colorbar;
    
    xlabel('Yaw Vel (deg/s)');
    ylabel('Fwd Vel (mm/s)');
    
    xlim([-200 200]);
    ylim([-5 15]);

end
