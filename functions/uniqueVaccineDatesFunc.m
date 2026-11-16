function [uniqueRows,isError] = uniqueVaccineDatesFunc(vaccineCoverageTable, dataFolder)
%set error flag false
isError = false;

% Find unique report dates in a stable order and get the index of their first occurrence
[uniqueElements, firstIdx, ~] = unique(vaccineCoverageTable, 'stable');
uniqueRows = vaccineCoverageTable;
if height(uniqueElements) ~= height(vaccineCoverageTable)

    % Construct filename to save the data in a specified folder
    filename = fullfile(getExperimentDir(), dataFolder, 'vaccineCoverage_Unique.txt');

    % Extract rows from the original table corresponding to the first occurrences of the unique elements
    uniqueRows = vaccineCoverageTable(firstIdx, :);

    % Save the unique rows to a text file
    writetable(uniqueRows, filename, 'Delimiter', 'tab');  % Using tab as a delimiter for better structure in text file

    % Set error falg to true
    isError = true;

    % Display success message and the unique elements
    disp('Error: the vaccine dates are not unique!')
    disp('Unique elements saved to file:');
    disp(filename);
end
end