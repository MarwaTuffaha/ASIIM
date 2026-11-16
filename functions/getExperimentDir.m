function fullExperimentDir = getExperimentDir(baseDir)
try
    load(fullfile(baseDir,'config.mat'),...
        'experimentDir');
catch ME
    error('Failed to load required variables: %s', ME.message);
end
fullExperimentDir = experimentDir;
end
