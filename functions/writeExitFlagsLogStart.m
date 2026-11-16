function exitFlagsFileID = writeExitFlagsLogStart(s)
    % Create or open the text file in the specified folder to write exit flags

    fullPathToFile = fullfile(s.folder, s.diagnosticFolder,  'exit_flags.txt');

    % Check if the directory exists, if not, create it
    [logPath, ~, ~] = fileparts(fullPathToFile);
    if ~exist(logPath, 'dir')
        mkdir(logPath);  % Create the directory if it does not exist
    end

    % Attempt to open the file
    exitFlagsFileID = fopen(fullPathToFile, 'a');
    
    % Check if the file was successfully opened
    if exitFlagsFileID == -1
        error('Failed to open log file: %s', fullPathToFile);
    end

    exitFlagsFileID = fopen(fullPathToFile, 'a'); % 'a' for append mode
    
    % Write the diagnostic log start message
    fprintf(exitFlagsFileID, '--------------------\n');
    fprintf(exitFlagsFileID, '### Diagnostic log started on %s ###\n', ...
        string(datetime('now')));
    fprintf(exitFlagsFileID, 'index     windowStart     windowEnd       exitFlag\n');
end
