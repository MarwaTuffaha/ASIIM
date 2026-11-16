function loopList = determineLoopList(s, fp, mitList)
    % Handle the window loop end
    if fp.windowLoopEnd >= mitList(end)
        fp.windowLoopEnd = mitList(end - 1);
    end

    % Find the start and end indices for mitagation
    mitIndexStart = find(fp.windowLoopStart <= mitList, 1);
    mitIndexEnd = find(fp.windowLoopEnd < mitList, 1);

    % Create the loop list
    loopList = mitList(mitIndexStart:fp.windowStep:mitIndexEnd);

    if s.hospital ==1
        dates=s.hospDates;
    else
        dates=s.incidenceDates;
    end

    % Correct the loop list if it overlaps with hospDates
    while dates(1) > s.dateNumList(loopList(2)) && length(loopList) > 3
        loopList = [loopList(1), loopList(3:end)];
    end

    % Check for date overlap error
    if dates(1) > s.dateNumList(loopList(2))
        disp('date overlap error');
    end

    % Check length of loopList
    if length(loopList) < 3
        disp('length of loopList <= 2 error');
    end
end