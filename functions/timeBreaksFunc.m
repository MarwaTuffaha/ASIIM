function [timeBreaks] = timeBreaksFunc(dateList)
%timeBreaksFunc Find the mitagation windows in simulation time. (i.e. from
%0.)
%   Detailed explanation goes here
day0 = datenum(datetime(...
    [dateList(1,3),dateList(1,2),dateList(1,1)]));
len = length(dateList);
timeBreaks = ones(1,len);
for index = 1:len
    timeBreaks(index) = datenum(datetime(...
        [dateList(index,3),dateList(index,2),dateList(index,1)]...
        ))-day0;
end
end


