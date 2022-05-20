% quick and dirty plot of single leg trials for given fly

% load file first

% some params
% numTrials = 5;
whichTrials = 16:20;
numParams = 2;
numDurs = length(durs);
% numDurs = 1;
legInd = 1:6;
whichND = 1.3;

whichNDInd = find(whichND == NDs);

% number of rows and columns of subplots
numRows = numTrials;
numCols = numParams * numDurs;

% plot parameters
yScale = [-1 1];


figure;
hold on;

counter = 1;

for i = whichTrials
    for k = 1:numDurs
        for j = 1:numParams
            subplot(numRows,numCols,counter);
            line([0 0],yScale,'Color','k');
            line([durs(k) durs(k)],yScale,'Color','k');
            hold on;

            if (j==1)
                thisParamVals = squeeze(legTrackOpto.srnfLegX.reps(whichNDInd,k,legInd));
                if (i==whichTrials(1))
                    title(sprintf('Leg X position, duration %.1f',durs(k)));
                end
            elseif (j==2)
                thisParamVals = squeeze(legTrackOpto.srnfLegY.reps(whichNDInd,k,legInd));
                if (i==whichTrials(1))
                    title(sprintf('Leg Y position, duration %.1f',durs(k)));
                end
            end
            for l = 1:length(legInd)
                thisLegVals = thisParamVals{legInd(l)};
                plot(legTrackOpto.durTs{k},thisLegVals(i,:));
                ylim(yScale);
            end
            
            counter = counter + 1;
        end
    end
end