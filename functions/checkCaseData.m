function errors = checkCaseData(data, obj)
% checkCaseData Checks for issues in the provided data and logs errors
%
% Syntax:
%   errors = checkCaseData(data, obj)
%
% Description:
%   This function checks if the provided data is non-empty and contains no NaN values.
%   If any issues are found, it writes the errors to a log file in the specified directory
%   and sets an error flag.
%
% Input Arguments:
%   data - The data to be checked (numeric array)
%   obj  - An object containing the error log file path
%
% Output Arguments:
%   errors - A flag indicating if any errors were found (0 if no errors, 1 if errors)
%
% Detailed Explanation:
%   1. The function defines the path for the log file based on the object's properties.
%   2. It initializes an error flag to 0 (no errors).
%   3. The function checks if the data is non-empty and contains no NaN values.
%      If any issues are found, it sets the error flag to 1, writes the error to the log file,
%      and displays an error message.
%
% Example:
%   errors = checkCaseData(hospitalData, simulationObject)
%
% Authors:
%   Your Name (your.email@example.com)
%   Date: July 12, 2024

    % Define the path for the log file
    logFilePath = obj.errorLogFilePath;
    errors = 0;  % Initialize the error flag

    try
        % Check if data is non-empty
        if isempty(data)
            errors = 1;  % Set error flag
            writeError(logFilePath, 'The provided data is empty.');
            error('The provided data is empty.');
        end

        % Check for NaN values in the data
        if any(isnan(data), 'all')
            errors = 1;  % Set error flag
            writeError(logFilePath, 'The provided data contains NaN values.');
            error('The provided data contains NaN values.');
        end
    catch e
        % Display error message in the console
        disp(['Error: ', e.message]);
    end
end