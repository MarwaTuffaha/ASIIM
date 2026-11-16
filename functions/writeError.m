function writeError(filePath, errMsg)
    % Open the file for appending and write the error message
    fid = fopen(filePath, 'a');
    if fid == -1
        error('Unable to open the log file.');
    end
    fprintf(fid, '%s - %s\n', datetime('now'), errMsg);
    fclose(fid);
end