function [contactsData,contactsDates] = contactsIn(path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fullPath = fullfile(path);
contactRaw = importdata(fullPath);

regions = unique(contactRaw.textdata(2:end,1));

len = length(regions);
for i = 1:len
%contacts(regexprep(provListUnique{i}, ' ', '')).region = regions{i}
contactsDates.(regexprep(regions{i}, ' ', '')) = datenum(datetime(...
    contactRaw.data(strcmp(contactRaw.textdata(2:end,1),regions{i}),1)...
    ,'ConvertFrom','yyyymmdd'));
contactsData.(regexprep(regions{i}, ' ', '')) =...
    contactRaw.data(strcmp(contactRaw.textdata(2:end,1),regions{i}),2:4);
end



%regexprep(provListUnique{i}, ' ', '')


end

