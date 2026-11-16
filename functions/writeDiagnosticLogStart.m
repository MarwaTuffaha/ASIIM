function diagnosticFileID = writeDiagnosticLogStart(s, fp)
    % Write the starting message of the diagnostic log.
    %
    % Parameters:
    %   - s: struct containing folder and diagnosticFolder
    %   - fp: struct containing windowsize and timeLimit

    % Get the number of available cores
    numCores = feature('numCores');
    
    % MATLAB version
    matlabVersion = version;

    % Operating system information
    osInfo = computer;

    % Check folder exists
    fullPath = fullfile(s.experimentDir, s.folder, s.diagnosticFolder);

    if ~exist(fullPath, 'dir')
        % Folder does not exist, create it
        mkdir(fullPath);
        disp(['Created directory: ', fullPath]);
    end

    % Create or open the text file in the specified folder to write diagnostics
    fullPathToFile = fullfile(s.experimentDir,s.folder,...
        s.diagnosticFolder, 'diagnostic_log.txt');

    diagnosticFileID = fopen(fullPathToFile, 'a'); % 'a' for append mode

    % Write the diagnostic log start message
    fprintf(diagnosticFileID,'--------------------\n');
    fprintf(diagnosticFileID, '### Diagnostic log started on %s ###\n', ...
        string(datetime('now')));
    fprintf(diagnosticFileID, 'MATLAB Version: %s\n', matlabVersion);
    fprintf(diagnosticFileID, 'Operating System: %s\n', osInfo);
    fprintf(diagnosticFileID, 'Number of cores available for GA: %d\n', numCores);
    fprintf(diagnosticFileID,'Window size: %d\n', fp.windowSize);
    fprintf(diagnosticFileID,'Window step: %d\n', fp.windowStep);
    fprintf(diagnosticFileID,'Time limit per step: %d\n', fp.timeLimit);
    fprintf(diagnosticFileID, 'Overall fitting window from: %s (%d) to: %s (%d)\n', ...
        datestr(s.dateNumList(fp.windowLoopStart,:)), fp.windowLoopStart, ...
        datestr(s.dateNumList(fp.windowLoopEnd,:)), fp.windowLoopEnd);
    fprintf(diagnosticFileID,'Curvature Penalty Percent: %f\n', fp.curvaturePenaltyPercent);
    fprintf(diagnosticFileID,'--------------------\n');
end
