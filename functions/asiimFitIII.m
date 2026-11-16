function asiimFitIII(s)
% asiimFitII Performs rolling-window GA fitting of model parameters
%            to incidence, hospital, or both, depending on s.hospital.
%
% Syntax:
%   asiimFitII(s)
%
% Input Arguments:
%   s  - An object containing simulation data/model and configuration
%
% Description:
%   s.hospital determines the fitting mode:
%       0 -> incidence only  (use s.fitProportions_inc, e.g. [0.6, 1])
%       1 -> hospital only   (use s.fitProportions_hosp, e.g. 0.2)
%       2 -> both            (fit both incidence & hospital simultaneously)
%
%   The rest of the code performs a rolling-window fit using a genetic
%   algorithm, logging diagnostics, applying a curvature penalty, etc.

    %% 1) Determine the fitting windows
    mitList   = pickMitagationWindows(s, 0, 0);
    loopList  = determineLoopList(s, s, mitList);

    %% 2) Pre-fit diagnostic plot
    s.simulate();
    s.simulate([loopList(1), loopList(end)]);
    s.fitPlot(true, [loopList(1), loopList(end)], true);

    %% 3) Open diagnostic logs
    diagnosticFileID = writeDiagnosticLogStart(s, s);
    exitFlagsFileID  = writeExitFlagsLogStart(s);

    %% 4) Save original k-values
    originalkVal = s.kval;

    %% fit proportions
    fitProp.fitProportions_hosp = s.fitProportions_hosp;
    fitProp.fitProportions_inc = s.fitProportions_inc;

    %% 5) Prepare the data based on s.hospital
    %    We'll store them in separate arrays. Some may remain empty if unused.
    incDates = s.incidenceDates;
    incCases = s.incidenceCases;
    hospDates = s.hospDates;
    hospCases = s.hospCases;

    lag = s.hospitalReportingLag;

    incDateNS  = [];
    incDataNS  = [];
    hospDateNS = [];
    hospDataNS = [];

    switch s.hospital
        case 0  % Incidence-only
            fprintf(diagnosticFileID, 'Fitting mode: Incidence ONLY (s.hospital=0)\n');

            incDateNS = incDates(incDates < s.stopDate);
            incDataNS = [incCases(incDates < s.stopDate); incCases(end) ];

        case 1  % Hospital-only
            fprintf(diagnosticFileID, 'Fitting mode: Hospital ONLY (s.hospital=1)\n');

            % Keep observed hospital data on its reported calendar dates.
            % Apply hospitalReportingLag to the model prediction instead.
            hospDateNS = hospDates(hospDates < s.stopDate);
            hospDataNS = hospCases(hospDates < s.stopDate, :);

        case 2  % Both
            fprintf(diagnosticFileID, 'Fitting mode: BOTH Incidence & Hospital (s.hospital=2)\n');

            incDateNS = incDates(incDates < s.stopDate);
            incDataNS = [ incCases(incDates < s.stopDate);  incCases(end) ];

            % Keep observed hospital data on its reported calendar dates.
            % Apply hospitalReportingLag to the model prediction instead.
            hospDateNS = hospDates(hospDates < s.stopDate);
            hospDataNS = hospCases(hospDates < s.stopDate, :);

        otherwise
            fprintf(diagnosticFileID, 'ERROR: s.hospital = %d is invalid.\n', s.hospital);
            error('Unrecognized s.hospital value = %d', s.hospital);
    end

    % Quick data consistency check
    diagnosticMsg = checkDataConsistency(s);
    fprintf(diagnosticFileID, '%s\n', diagnosticMsg);

    % Start timing
    thisFitTiming = 0;

    %% 6) Main rolling-window loop
    len = length(loopList);
    lastWindowIndex = 1;
    for i = 1:len
        if loopList(i) < (loopList(end) - s.windowSize)
            lastWindowIndex = i;
        end
    end
    % Halt if loop size is insufficient for fitting
    if lastWindowIndex < 2
        fprintf(diagnosticFileID, 'Fit failed to start: Increase windowLoopEnd-windowLoopStart\n');
        error('Fit failed to start, insufficient loop size.');
    end

    for index = 1:lastWindowIndex
        % Log the current fit loop index
        fprintf(diagnosticFileID, 'Fit loop index: %d\n', index);
        tic;
        
        % Determine the start and end indices for the fitting window
        indexStart = loopList(index);
        indexEnd   = loopList(index) + s.windowSize;

        % Log the fitting window details
        fprintf(diagnosticFileID, ...
            'Fitting window from: %s (%d) to: %s (%d)\n', ...
            datestr(s.dateNumList(indexStart,:)), indexStart, ...
            datestr(s.dateNumList(indexEnd,:)),   indexEnd);

        %----------------------------------------------------------------------
        % Curvature penalty initialization
        %----------------------------------------------------------------------
        if index == 1
            kvalStart = s.kval;
        else
            kvalStart = kvalMin;
        end

        % We’ll get the “raw” curvature penalty from user
        % curvature penalty (i.e. Regularization Hyperparameter)
        %With the assumption that the public changes in the level of mitagation
        %is 'sticky' in the economic sense. Rapid changes to the fit k-vlaue are
        %penalized using an approximation of curvature of the k-value vector.
        %Choose a value that is small relative to the error.
        %NOTE: this is also set after obj is normalized.
        curvaturePenalty = s.curvaturePenalty;

        %----------------------------------------------------------------------
        % Identify sub-vector of k-values for GA
        %----------------------------------------------------------------------
        mitIndexStart = find(indexStart <= mitList, 1);
        mitIndexEnd   = find(indexEnd   >= mitList, 1, "last");

        initKval = kvalStart(mitList(mitIndexStart:mitIndexEnd));
        % e.g., omit the last element if your logic says so
        %initKval = initKval(1:end-1);
        lenKvalStart = length(initKval);

        %----------------------------------------------------------------------
        % Create initial GA population
        %----------------------------------------------------------------------
        initPopSize = s.initPopSize;
        initPop = zeros(initPopSize, lenKvalStart);
        initPop(1,:) = initKval;

        for ii = 1:lenKvalStart
            initPop(2:end, ii) = normrnd(initKval(ii), 0.8,1,initPopSize-1);
        end
        % Clamp out-of-bounds to [s.lowerBound, s.upperBound]
        initPop(initPop < s.lowerBound) = 1;
        initPop(initPop > s.upperBound) = 1;

        %----------------------------------------------------------------------
        % Define objective function: first, evaluate at initKval to scale penalty
        %----------------------------------------------------------------------
        tempObj = @(x) forceFitObjFastFunc2_smooth( ...
            x, s, ...
            incDateNS, incDataNS, ...
            hospDateNS, hospDataNS, ...
            indexStart, indexEnd, ...
            mitIndexStart, mitIndexEnd, mitList, ...
            false, ...          % monitorPlot = false for initial call
            curvaturePenalty, ... 
            1, ...              % objStart placeholder
            s.fitProportions_inc, ...  % e.g. [0 0 0.6 1]
            s.fitProportions_hosp ...  % e.g. [0 0 0 0.2]
        );

        objStart = tempObj(initKval);
        % Rescale curvature penalty based on initial SSE
        curvaturePenalty = objStart * s.curvaturePenaltyPercent;

        % Now define the real objective handle
        f = @(x) forceFitObjFastFunc2_smooth( ...
            x, s, ...
            incDateNS, incDataNS, ...
            hospDateNS, hospDataNS, ...
            indexStart, indexEnd, ...
            mitIndexStart, mitIndexEnd, mitList, ...
            s.monitorPlot, curvaturePenalty, objStart, ...
            s.fitProportions_inc, s.fitProportions_hosp);

        fprintf(diagnosticFileID, ...
            'Objective Value at Starting Point: %f\n', objStart);
        fprintf(diagnosticFileID, ...
            'Curvature Penalty: %f\n', curvaturePenalty);

        %----------------------------------------------------------------------
        % GA configuration
        %----------------------------------------------------------------------
        populationSize = lenKvalStart;
        lb = s.lowerBound * ones(1, lenKvalStart);
        ub = s.upperBound * ones(1, lenKvalStart);

        A = [];
        b = [];
        Aeq = [];
        beq = [];

        % % Example GA Output function (optional)
        % nameStruct.region   = s.region;
        % nameStruct.scenario = s.simNum;
        % gaOutNest = @(options, state, flag) gaOutputPopulation(options, state, flag, nameStruct);

        options = optimoptions('ga', ...
            'UseParallel', s.useParallel, ...
            'InitialPopulationMatrix', initPop, ...
            'PlotFcn', {@gaplotbestf, @gaplotdistance, @gaplotbestindiv, @gaplotselection}, ...
            'Display', 'iter', ...
            'MaxStallTime', s.maxStallTime, ...
            'MaxStallGenerations', 20, ...
            'PopulationSize', s.populationFactor * populationSize, ...
            'MaxGenerations', s.maxGenerations, ...
            'MaxTime', s.timeLimit * 60 ...
            ... % 'OutputFcn', gaOutNest,  % Uncomment if you want to track the GA population
        );

        %----------------------------------------------------------------------
        % Run GA
        %----------------------------------------------------------------------
        [kvalFMinShort, fval, exitflag, ~] = ga(f, lenKvalStart, A, b, Aeq, beq, lb, ub, [], options);

        fprintf(diagnosticFileID, 'GA exit flag: %d\n', exitflag);
        fprintf(diagnosticFileID, 'Best fit objective function value: %f\n', fval);

        % record exit flags
        fprintf(exitFlagsFileID, '%d     %s     %s       %d\n', ...
            index, datestr(s.dateNumList(indexStart,:)), ...
            datestr(s.dateNumList(indexEnd,:)), exitflag);

        %----------------------------------------------------------------------
        % Expand k-values back to full s.kval
        %----------------------------------------------------------------------
        kvalFull = kvalStart;
        kvalMin  = kValueExpand(mitIndexStart, mitIndexEnd, mitList, kvalFMinShort, kvalFull);

        % Update s.kval
        s.kval = kvalMin;
        writematrix(kvalMin, s.kValueHandle);

        %----------------------------------------------------------------------
        % Diagnostic plot
        %----------------------------------------------------------------------
        s.simulate();
        s.simulate([indexStart, indexEnd]);
        createDiagnosticPlot(s, originalkVal, kvalMin, indexStart, indexEnd, index);
        close
        disp(loopList(index))

        %----------------------------------------------------------------------
        % Timing and break on NaN
        %----------------------------------------------------------------------
        fitTimming     = toc;
        thisFitTiming  = thisFitTiming + fitTimming;
        s.fitTime      = s.fitTime + fitTimming;

        fprintf(diagnosticFileID, 'Time for this loop: %f seconds\n', fitTimming);

        if isnan(fval)
            fprintf(diagnosticFileID, ...
                'Exiting loop. Objective function returned NaN.\n');
            break;
        end
    end

    %% 7) Final simulation & plot
    s.simulate();
    s.fitPlot();
    plotBothData(s)

    %% 8) Close logs
    writeDiagnosticLogEnd(diagnosticFileID, f, initKval, objStart, ...
                          thisFitTiming, fitTimming);
    fclose(exitFlagsFileID);
end
        
