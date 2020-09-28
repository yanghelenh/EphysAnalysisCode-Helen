% loadMetadataSpreadsheet.m
%
% Function to load information in metadata spreadsheet into a data struct, 
%  to allow selection of specific subset of data for analysis.
%
% INPUT:
%   none - but brings up GUI to allow user to select metadata spreadsheet
%
% OUTPUT:
%   metaDat - struct of metadata
% 
% CREATED: 9/27/20 - HHY
%
% UPDATED: 9/27/20 - HHY
%

function metaDat = loadMetadataSpreadsheet()

    % get metadata spreadsheet path
    disp('Select metadata spreadsheet for experiment');
    [sprdshtName, sprdshtPath] = uigetfile('*.xlsx', ...
        'Select metadata spreadsheet');
    sprdshtFullPath = [sprdshtPath filesep sprdshtName];
    
    % get information about spreadsheet
    opts = detectImportOptions(sprdshtFullPath);
    
    % change options to import relevant metadata (exclude freeform comments
    %  fields) in correct format
    opts.SelectedVariableNames = opts.VariableNames([1:11, 13:20]);
    
    opts = setvaropts(opts,'Exclude','FillValue',0);
    
    % read in spreadsheet as table
    sprdsht = readtable(sprdshtFullPath, opts);
    
    % convert to struct
    for i = 1:size(sprdsht,2)
        metaDat.(sprdsht.Properties.VariableNames{i}) = table2cell(...
            sprdsht(:,i));
    end

    % convert FlyID, Age, Exclude to normal arrays
    metaDat.Age = cell2mat(metaDat.Age);
    metaDat.FlyID = cell2mat(metaDat.FlyID);
    metaDat.Exclude = cell2mat(metaDat.Exclude);

end
