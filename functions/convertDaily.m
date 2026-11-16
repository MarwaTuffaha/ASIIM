function infectedDaily = convertDaily(rp)
% convertDaily Computes the daily number of newly infected individuals
%
% Syntax:
%   infectedDaily = convertDaily(rp)
%
% Description:
%   This function calculates the daily count of new infections from cumulative 
%   infection data. It processes the data for each compartment and age group, 
%   returning the daily infections as the difference between consecutive days.
%
% Input Arguments:
%   rp - A structure array containing the cumulative infection data across 
%        different compartments and age groups.
%
% Output Arguments:
%   infectedDaily - A 3D array representing the daily new infections. The 
%                   dimensions correspond to days, age groups, and compartments.

len = length(rp(4).C);

infectedDaily = zeros(len-1,16,8);%%1+5,2+6,3+7,4+8 are combined in C % Not anymore!
for k = [2:4, 6:8]
    for j = 1:16
        infectedDaily(1:end,j,k) = rp(k).C(2:end,j)-rp(k).C(1:end-1,j);
    end
end
end