function [objScaled] = forceFitObjFastFunc2_smooth(...
    kvalIn, s, ...
    incDateNS, incDataNS, ...
    hospDateNS, hospDataNS, ...
    indexStart, indexEnd, ...
    mitIndexStart, mitIndexEnd, mitList, ...
    plotBool, curvaturePenalty, objStart, ...
    fitPropInc, fitPropHosp)
% forceFitObjFastFunc2
%
% Objective function for asiimFitIII, handling:
%   s.hospital = 0 -> Incidence only
%   s.hospital = 1 -> Hospital only
%   s.hospital = 2 -> Both incidence & hospital
%
% This version:
%   1) Runs the model with the provided k-values.
%   2) Converts daily or keeps the data as needed.
%   3) Uses proportionToFitFunc(...) to get incModel or hospModel (2D).
%   4) Aggregates (sums) across the 2nd dimension to yield a 1D time series.
%   5) Compares to 1D observed data, computing a single RMSE-like metric.
%   6) Adds a curvature penalty if k-values change rapidly.
%   7) Returns a SINGLE final scalar = (objInc + objHosp + penalty) / objStart
%
% If incModel or hospModel is more than 2D, or if you need a different
% aggregation method (e.g., weighted sums), adjust the "sum(...,2)" line.

    %% 1) Initialize the objective
    obj = 0;

    %% 2) Expand GA/pattern-search sub-vector (kvalIn) into s.kval
    oldKval = s.kval;
    s.kval  = kValueExpand(mitIndexStart, mitIndexEnd, mitList, kvalIn, oldKval);

    %% 3) Curvature penalty (second derivative of k-values)
    if length(kvalIn) >= 3
        ddk      = diff(kvalIn, 2).^2;         % second derivative squared
        curvTerm = sum(ddk) / length(kvalIn);  % average or sum
        obj      = obj + (curvaturePenalty * curvTerm);
    end

    %% 4) Run the model/solver with these parameters
    %s.simulate
    [time, yout] = s.fitSimulate([mitList(mitIndexStart), mitList(mitIndexEnd)]);
    
    % IMPORTANT:
    % getPrevLognormal/getPrevPERT call getPrediction(s), and getPrediction(s)
    % uses s.yout and s.time. Therefore these fields must be updated here.
    s.time = time;
    s.yout = yout;
    
    %% 5) Post-process
    rp = postProcess(yout);

    %% 6) Convert to daily if needed
    %    (If hospital data is also daily or if you have different flags, adapt)
    % if s.hospital == 1
    %     dataTypeLabel = s.dataTypeLabel_hosp;  % hospital-only
    % else
    %     dataTypeLabel = s.dataTypeLabel_inc;   % incidence or both
    % end

    if or(or(s.hospital == 0 && s.dataTypeLabel_inc==0,s.hospital == 1 && s.dataTypeLabel_hosp==0),and(s.hospital == 2 && s.dataTypeLabel_inc==0, or(s.dataTypeLabel_hosp==0, s.dataTypeLabel_hosp==2))) % data is daily incidence and/or daily hospital admissions
        % Convert to daily new cases/hospitalizations
        modelSolution = convertDaily(rp);
        % time is typically one less if we do daily differences
        time = time(1:end-1);

        % If fitting incidence together with hospital prevalence, also build
        % the LOS-based hospital prevalence prediction for the hospital part.
        if s.hospital == 2 && s.dataTypeLabel_hosp==2
            if s.LOS_dist=="lognormal"
                [t_hosp_prev, midSol_hospPrev] = getPrevLognormal(s, s.LOS_Q1, s.LOS_median, s.LOS_Q3);
            elseif s.LOS_dist=="PERT"
                [t_hosp_prev, midSol_hospPrev] = getPrevPERT(s, s.LOS_Q1, s.LOS_median, s.LOS_Q3);
            else
                error("Unknown LOS distribution: %s", s.LOS_dist);
            end
        end

    elseif s.hospital == 0 && s.dataTypeLabel_inc==2  % We only have incidence prevalence data     
        modelSolution = zeros(size(rp(1).I,1), size(rp(1).I,2), 8);
        for i = 1:8
            modelSolution(:,:,i) = rp(i).I;
        end
        modelSolution = modelSolution(1:end,:,:);

    elseif s.hospital == 1 && s.dataTypeLabel_hosp==2  % We only have hospital prevalence data

        if s.LOS_dist=="lognormal"
            [t_hosp_prev, midSol_hospPrev] = getPrevLognormal(s, s.LOS_Q1, s.LOS_median, s.LOS_Q3);
        elseif s.LOS_dist=="PERT"
            [t_hosp_prev, midSol_hospPrev] = getPrevPERT(s, s.LOS_Q1, s.LOS_median, s.LOS_Q3);
        else
            error("Unknown LOS distribution: %s", s.LOS_dist);
        end

    else
        error("Code does not deal with these types/combinations of data. Either incidence and hospital daily, or only incidence prevalence, or only hospital prevalence")
    end

    %% 7) Identify overlapping time points
    day0Simulation = s.dateNumList(1);

    modelIncDates  = day0Simulation + time;
    modelHospDates = day0Simulation + time + s.hospitalReportingLag;

    bothIncData = [];
    bothSimInc  = [];
    if ~isempty(incDateNS)
        bothIncData = ismember(incDateNS,  modelIncDates);
        bothSimInc  = ismember(modelIncDates, incDateNS(bothIncData));
    end

    bothHospData = [];
    bothSimHosp  = [];
    if ~isempty(hospDateNS)
        bothHospData = ismember(hospDateNS, modelHospDates);
        bothSimHosp  = ismember(modelHospDates, hospDateNS(bothHospData));
    end

    %% 8) Initialize partial objectives (for incidence & hospital)
    objInc  = 0;
    objHosp = 0;

    %% 9) Depending on s.hospital, compute single-scalar partial errors
    switch s.hospital

        case 0
            % ======= Incidence Only =======
            if any(bothSimInc)

                midDaily = proportionToFitFunc(fitPropInc, modelSolution, bothSimInc);
                midSol_inc = (midDaily(bothSimInc));
                [~, smoothData_inc] = sumAndSmoothSubset(incDataNS, bothIncData, 7);
                smoothData_inc = abs(smoothData_inc); % avoid negative values

                residInc = zeros(length(smoothData_inc),1);
                for i = 1:length(smoothData_inc)
                    if smoothData_inc(i) ~= 0 || midSol_inc(i) ~= 0
                        %residInc(i)=(smoothData_inc(i)-midSol_inc(i))./max(eps,midSol_inc(i)); Relative error
                        residInc(i) = (smoothData_inc(i)-midSol_inc(i))./max(eps,(smoothData_inc(i)+midSol_inc(i))/2); % percent error
                    end
                end
            
                objInc = sum(abs(residInc))/max(1,sum(bothSimInc));

            end

        case 1
            % ======= Hospital Only =======
            % s.hospital == 1 means the objective is based only on hospital data.
            % s.dataTypeLabel_hosp decides whether this hospital data is daily
            % admissions (0) or hospital prevalence (2).

            if s.dataTypeLabel_hosp==0 % data is daily hosp

                if any(bothSimHosp)

                    midDaily = proportionToFitFunc(fitPropHosp, modelSolution, bothSimHosp);
                    midSol_hosp = midDaily(bothSimHosp);

                    if s.hospCasesSmoothLabel==0
                        [~, smoothData_hosp] = sumAndSmoothSubset(hospDataNS, bothHospData, 7);
                    elseif s.hospCasesSmoothLabel==1 % data is already smoothed
                        smoothData_hosp = sumColumnsIfNeeded(hospDataNS(bothHospData, :));
                    else
                        error('Unrecognized s.hospCasesSmoothLabel=%d', s.hospCasesSmoothLabel);
                    end

                else
                    error('No overlap between observed daily hospital dates and model dates.');
                end

            elseif s.dataTypeLabel_hosp==2 % data is hosp prevalence

                % Hospital prevalence from the LOS function uses its own date vector.
                % Therefore, align observed hospital prevalence to t_hosp_prev,
                % not to day0Simulation + time.
                midDaily = midSol_hospPrev;
                hospPrevDates = t_hosp_prev;

                % Use rounded dates to avoid tiny datenum precision mismatches.
                hospDateNS_round    = round(hospDateNS);
                hospPrevDates_round = round(hospPrevDates);

                bothHospPrevData = ismember(hospDateNS_round, hospPrevDates_round);
                bothHospPrevSim  = ismember(hospPrevDates_round, hospDateNS_round(bothHospPrevData));

                if ~any(bothHospPrevData) || ~any(bothHospPrevSim)
                    disp('--- Hospital prevalence fitting date debug ---')
                    disp(['min(hospDateNS) = ', num2str(min(hospDateNS))])
                    disp(['max(hospDateNS) = ', num2str(max(hospDateNS))])
                    disp(['min(hospPrevDates) = ', num2str(min(hospPrevDates))])
                    disp(['max(hospPrevDates) = ', num2str(max(hospPrevDates))])
                    disp(['min(t_hosp_prev) = ', num2str(min(t_hosp_prev))])
                    disp(['max(t_hosp_prev) = ', num2str(max(t_hosp_prev))])
                    disp(['day0Simulation = ', num2str(day0Simulation)])

                    error(['No overlap between observed hospital prevalence dates and ' ...
                           'LOS-based model prevalence dates. Check t_hosp_prev and hospDateNS.']);
                end

                midSol_hosp = midDaily(bothHospPrevSim);

                if s.hospCasesSmoothLabel==0
                    [~, smoothData_hosp] = sumAndSmoothSubset(hospDataNS, bothHospPrevData, 7);
                elseif s.hospCasesSmoothLabel==1 % data is already smoothed
                    smoothData_hosp = sumColumnsIfNeeded(hospDataNS(bothHospPrevData, :));
                else
                    error('Unrecognized s.hospCasesSmoothLabel=%d', s.hospCasesSmoothLabel);
                end

            else
                error('Hospital-only calibration requires s.dataTypeLabel_hosp to be 0 daily or 2 prevalence.');
            end

            smoothData_hosp = abs(smoothData_hosp); % avoid negative values

            % Force column vectors
            smoothData_hosp = smoothData_hosp(:);
            midSol_hosp = midSol_hosp(:);

            % Strict alignment check. If this fails, the model and data are
            % not being compared on the same dates.
            if length(smoothData_hosp) ~= length(midSol_hosp)
                error('Hospital alignment problem: data length = %d, model length = %d', ...
                    length(smoothData_hosp), length(midSol_hosp));
            end

            residHosp = zeros(length(smoothData_hosp),1);
            for i = 1:length(smoothData_hosp)
                if smoothData_hosp(i) ~= 0 || midSol_hosp(i) ~= 0
                    %residHosp(i)=(smoothData_hosp(i)-midSol_hosp(i))./max(eps,midSol_hosp(i));
                    residHosp(i) = abs(smoothData_hosp(i)-midSol_hosp(i))./max(eps,(smoothData_hosp(i)+midSol_hosp(i))/2);
                end
            end

            objHosp = sum(abs(residHosp))/max(1,length(residHosp));

        case 2
            % ======= Both Incidence & Hospital =======
            % 1) Incidence Part
            if any(bothSimInc)

                midDaily = proportionToFitFunc(fitPropInc, modelSolution, bothSimInc);
                midSol_inc = (midDaily(bothSimInc));
                [~, smoothData_inc] = sumAndSmoothSubset(incDataNS, bothIncData, 7);
                smoothData_inc = abs(smoothData_inc); % avoid negative values

                residInc = zeros(length(smoothData_inc),1);
                for i = 1:length(smoothData_inc)
                    if smoothData_inc(i) ~= 0 || midSol_inc(i) ~= 0
                        %residInc(i)=(smoothData_inc(i)-midSol_inc(i))./max(eps,midSol_inc(i));
                        residInc(i) = (smoothData_inc(i)-midSol_inc(i))./max(eps,(smoothData_inc(i)+midSol_inc(i))/2);
                    end
                end
            
                objInc = sum(abs(residInc))/max(1,sum(bothSimInc));                

            end

            % 2) Hospital Part
            % If this fitting window has incidence overlap but no hospital
            % overlap, ignore the hospital term for this window instead of
            % stopping the whole calibration.
            fitHospThisWindow = false;

            if s.dataTypeLabel_hosp==0 % data is daily hospital admissions

                if any(bothSimHosp)
                    fitHospThisWindow = true;
                    midDaily = proportionToFitFunc(fitPropHosp, modelSolution, bothSimHosp);
                    midSol_hosp = midDaily(bothSimHosp);

                    if s.hospCasesSmoothLabel==0
                        [~, smoothData_hosp] = sumAndSmoothSubset(hospDataNS, bothHospData, 7);
                    elseif s.hospCasesSmoothLabel==1 % data is already smoothed
                        smoothData_hosp = sumColumnsIfNeeded(hospDataNS(bothHospData, :));
                    else
                        error('Unrecognized s.hospCasesSmoothLabel=%d', s.hospCasesSmoothLabel);
                    end
                elseif ~any(bothSimInc)
                    error(['No overlap between observed case dates and model dates, ' ...
                           'and no overlap between observed daily hospital dates and lagged model hospital dates.']);
                end

            elseif s.dataTypeLabel_hosp==2 % data is hospital prevalence

                % Hospital prevalence uses the LOS-based date vector from
                % getPrevLognormal/getPrevPERT, which already inherits the
                % hospitalReportingLag through getPrediction.
                hospPrevDates = t_hosp_prev;

                hospDateNS_round    = round(hospDateNS);
                hospPrevDates_round = round(hospPrevDates);

                bothHospPrevData = ismember(hospDateNS_round, hospPrevDates_round);
                bothHospPrevSim  = ismember(hospPrevDates_round, hospDateNS_round(bothHospPrevData));

                if any(bothHospPrevData) && any(bothHospPrevSim)
                    fitHospThisWindow = true;
                    midDaily = midSol_hospPrev;
                    midSol_hosp = midDaily(bothHospPrevSim);

                    if s.hospCasesSmoothLabel==0
                        [~, smoothData_hosp] = sumAndSmoothSubset(hospDataNS, bothHospPrevData, 7);
                    elseif s.hospCasesSmoothLabel==1 % data is already smoothed
                        smoothData_hosp = sumColumnsIfNeeded(hospDataNS(bothHospPrevData, :));
                    else
                        error('Unrecognized s.hospCasesSmoothLabel=%d', s.hospCasesSmoothLabel);
                    end
                elseif ~any(bothSimInc)
                    error(['No overlap between observed case dates and model dates, ' ...
                           'and no overlap between observed hospital prevalence dates and LOS-based model prevalence dates.']);
                end

            else
                error('Simultaneous calibration requires s.dataTypeLabel_hosp to be 0 daily or 2 prevalence.');
            end

            if fitHospThisWindow
                smoothData_hosp = abs(smoothData_hosp);
                smoothData_hosp = smoothData_hosp(:);
                midSol_hosp = midSol_hosp(:);

                if length(smoothData_hosp) ~= length(midSol_hosp)
                    error('Hospital alignment problem: data length = %d, model length = %d', ...
                        length(smoothData_hosp), length(midSol_hosp));
                end

                residHosp = zeros(length(smoothData_hosp),1);
                for i = 1:length(smoothData_hosp)
                    if smoothData_hosp(i) ~= 0 || midSol_hosp(i) ~= 0
                        residHosp(i) = abs(smoothData_hosp(i)-midSol_hosp(i))./max(eps,(smoothData_hosp(i)+midSol_hosp(i))/2);
                    end
                end
            
                objHosp = sum(abs(residHosp))/max(1,length(residHosp));
            end                

        otherwise
            error('Unrecognized s.hospital=%d', s.hospital);
    end

    %% 10) Combine partial objectives + curvature penalty
    obj = obj + objInc + objHosp;

    %% 11) Scale final objective
    objScaled = obj * (1 / objStart);

end


function [sumData, smoothData] = sumAndSmoothSubset(dataMatrix, bothDataMask, smoothWindow)
% sumAndSmoothSubset
%   1) Subset dataMatrix by bothDataMask.
%   2) Sum columns if multi-column.
%   3) Smooth with a given window using LOESS.
%
    subsetData = dataMatrix(bothDataMask, :);
    sumData    = sumColumnsIfNeeded(subsetData);
    smoothData = smooth(sumData, smoothWindow, 'loess');
end


function mergedData = sumColumnsIfNeeded(dataIn)
% sumColumnsIfNeeded
%   If dataIn has multiple columns, sum across them.
%   Otherwise, return dataIn as is.
%
    if size(dataIn, 2) > 1
        mergedData = sum(dataIn, 2);
    else
        mergedData = dataIn;
    end
end