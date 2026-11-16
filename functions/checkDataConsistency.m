function diagnosticMsg = checkDataConsistency(s)
% checkDataConsistency Verifies the integrity and consistency of hospital data
%
% Syntax:
%   diagnosticMsg = checkDataConsistency(s)
%
% Description:
%   This function checks the consistency of the hospital data in the simulation structure `s`. 
%   It verifies that the date-related fields (`hospDates`, `stopDate`, `hospitalReportingLag`) 
%   are numeric, checks if `hospDates` are sorted, and ensures that the lengths 
%   of `hospCases` and `hospDates` are equal. The function returns a message indicating 
%   the results of these checks.
%
% Input Arguments:
%   s - A structure containing simulation data and configuration, including fields 
%       such as `hospDates`, `stopDate`, `hospitalReportingLag`, and `hospCases`.
%
% Output Arguments:
%   diagnosticMsg - A string containing the diagnostic message. It lists any detected 
%                   errors or confirms that all checks passed.
%
% Examples:
%   % Example 1: Checking data consistency for a simulation object
%   simulationData = load('simulationData.mat');
%   diagnosticMsg = checkDataConsistency(simulationData);
%   disp(diagnosticMsg); % Output the result of the data consistency check
%
% Authors:
%   David Dick dwd416@gmail.com
    diagnosticMsg = '';

    % Check if s.hospDates, s.stopDate, and s.hospitalReportingLag are numeric
    if ~isnumeric(s.hospDates) || ~isnumeric(s.stopDate) || ~isnumeric(s.hospitalReportingLag)
        diagnosticMsg = [diagnosticMsg, 'Error: One or more date or lag variables are not numeric.'];
    end
    
    % Check if s.hospDates is sorted
    if ~issorted(s.hospDates)
        diagnosticMsg = [diagnosticMsg, 'Error: s.hospDates is not sorted.'];
    end
    
    % Check if s.hospCases and s.hospDates have the same length
    if length(s.hospCases) ~= length(s.hospDates)
        diagnosticMsg = [diagnosticMsg, 'Error: Length of s.hospCases and s.hospDates are not the same.'];
    end

    % If no errors, write a message indicating success
    if isempty(diagnosticMsg)
        diagnosticMsg = 'Data consistency check passed.';
    end
end
