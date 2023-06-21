% sortLegStepsByOpto_allGUI.m
%
% Function that categorizes all leg steps by whether they're during a
%  opto stim bout or not.
% Calls sortLegStepsByOpto() but operates on multiple pData files, selected
%  through GUI
%
% INPUTS:
%   durs - vector of all durations of stimulation to consider
%   optoTime - length 2 vector where 1st element is time after opto starts
%       to begin counting step as during opto (as time in sec relative to
%       opto start time) and 2nd element is time before opto ends to stop
%       counting step as during opto (as time in sec relative to opto end
%       time)
%   notOptoTime - time in sec after opto stim turns off to not include 
%       in not opto category
%
% OUTPUTS:
%   none, but saves vector of category labels and key back into same pData
%       file
%
% CREATED: 12/6/22 - HHY
%
% UPDATED:
%   12/6/22 - HHY
%
function sortLegStepsByOpto_allGUI(durs, optoTime, notOptoTime)

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

        % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];
        
        % function call
        sortLegStepsByOpto(durs, optoTime, notOptoTime, ...
            pDataFullPath);
    end
end
    