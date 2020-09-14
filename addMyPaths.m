% addMyPaths.m
% 
% function to add relevant ephys data analysis paths, including subfolders
%
% Users: change the following paths to match those on your local computer

function addMyPaths()    
    %% Animal Part Tracker code
    APTpath = '/Users/hyang/Documents/MATLAB/APT';
    addpath(genpath(APTpath));

    %% Analysis code repository (EphysAnalysisCode-Helen)
    analysisPath = '/Users/hyang/Documents/EphysAnalysisCode-Helen';
    addpath(genpath(analysisPath));

    %% Folder containing metadata spreadsheet 
    metadataPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen';
    addpath(genpath(metadataPath));
    
    %% Python functions, add to python system path
    a2libPath = '/Users/hyang/Documents/EphysAnalysisCode-Helen/a2lib';
    P = py.sys.path;
    if count(P,a2libPath) == 0
        insert(P,int32(0),a2libPath);
    end

    %% Folder containing processed data (*_pdata.mat files)
    addpath(pDataDir());

end