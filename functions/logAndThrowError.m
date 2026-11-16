function logAndThrowError(obj, errorMessage)
    % Open log file to append error message
    logFile = fopen(obj.errorLogFilePath, 'a');
    if logFile == -1
        error('Unable to open the log file.');
    end
    
    % Write the error message with a timestamp
    fprintf(logFile, '%s - Error encountered \n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    fprintf(logFile, '- Error: %s\n', errorMessage);
    
    % Close the log file
    fclose(logFile);
    
    % Throw error after logging
    error(errorMessage);
end