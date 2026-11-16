function [residIncSum, residHospSum]=plotBothData(obj)

    %% 1) Run postProcess on obj.yout
    rp = postProcess(obj.yout);
    time = obj.time;  % Time vector from your simulation

    %------------------------------------------
    %% 2) Prepare date & data arrays for incidence/hospital
    %------------------------------------------
    dateData_hosp = obj.hospDates - obj.hospitalReportingLag;
    data_hosp     = obj.hospCases;
    dateData_inc  = obj.incidenceDates;
    data_inc      = obj.incidenceCases;

    % Validate case data
    if ~isempty(dateData_inc)
        checkCaseData(data_inc, obj);
    end
    if ~isempty(dateData_hosp)
        checkCaseData(data_hosp, obj);
    end
    
    %------------------------------------------
    %% 3) Align simulation time with observed data
    %------------------------------------------
    date0          = obj.dateList(1,:);
    dateList       = datenum(obj.dateList(:,3), obj.dateList(:,2), obj.dateList(:,1));
    day0Simulation = datenum(date0(3), date0(2), date0(1));

    bothData_hosp  = ismember(dateData_hosp, dateList(1):dateList(end));
    bothData_inc   = ismember(dateData_inc,  dateList(1):dateList(end));

    dateData_inc_common=dateData_inc(bothData_inc);
    dateData_hosp_common=dateData_hosp(bothData_hosp);

    bothSim_hosp = ismember(time(1:end-1) + day0Simulation- obj.hospitalReportingLag, dateData_hosp_common);
    bothSim_inc  = ismember(time(1:end-1) + day0Simulation, dateData_inc_common);

    t_hosp         = time(bothSim_hosp) + day0Simulation- obj.hospitalReportingLag;
    t_inc          = time(bothSim_inc)  + day0Simulation;

    %------------------------------------------
    %% 4) Build full daily/prevalence solution from rp
    %------------------------------------------
    if obj.hospital == 1
        dataTypeLabel = obj.dataTypeLabel_hosp;  % hospital-only
    else
        dataTypeLabel = obj.dataTypeLabel_inc;   % incidence or both
    end    

    if dataTypeLabel==0 % data is given daily
        % Convert to daily new cases/hospitalizations
        modelSolution = convertDaily(rp);
        % time is typically one less if we do daily differences
        time = time(1:end-1);

    elseif dataTypeLabel==2 % data is given as prevalence

        % For hospital prevalence, use LOS-based prevalence estimate.
        % For incidence prevalence, keep the old rp(i).I approach.
        if obj.hospital==1 && obj.dataTypeLabel_hosp==2

            if obj.LOS_dist=="lognormal"
                [t_hosp_prev, midSol_hospPrev] = getPrevLognormal(obj, obj.LOS_Q1, obj.LOS_median, obj.LOS_Q3);
            elseif obj.LOS_dist=="PERT"
                [t_hosp_prev, midSol_hospPrev] = getPrevPERT(obj, obj.LOS_Q1, obj.LOS_median, obj.LOS_Q3);
            else
                error("Unknown LOS distribution: %s", obj.LOS_dist);
            end

            % IMPORTANT:
            % t_hosp_prev is already absolute datenum. Do NOT add day0Simulation.
            hospPrevDates = t_hosp_prev;

            dateData_hosp_round = round(dateData_hosp);
            hospPrevDates_round = round(hospPrevDates);

            bothData_hosp = ismember(dateData_hosp_round, hospPrevDates_round);
            bothSim_hosp  = ismember(hospPrevDates_round, dateData_hosp_round(bothData_hosp));

            dateData_hosp_common = dateData_hosp(bothData_hosp);
            t_hosp = hospPrevDates(bothSim_hosp);

            if isempty(t_hosp) || ~any(bothSim_hosp)
                error(['No overlapping dates found between observed hospital prevalence ' ...
                       'and LOS-based hospital prevalence model output. Check t_hosp_prev, ' ...
                       'dateData_hosp, and hospitalReportingLag.']);
            end

        else
            modelSolution = zeros(size(rp(1).I,1), size(rp(1).I,2), 8);
            for i = 1:8
                modelSolution(:,:,i) = rp(i).I;
            end
            modelSolution = modelSolution(1:end,:,:);
        end
    end    

    %------------------------------------------
    %% 5) Sum & smooth for PLOTTING subsets (inc/hosp) 
    %    so we can plot raw vs. smoothed vs. model
    %------------------------------------------
    [sumData_inc, smoothData_inc]   = sumAndSmoothSubset(data_inc,  bothData_inc,  7);
    [sumData_hosp, smoothData_hosp] = sumAndSmoothSubset(data_hosp, bothData_hosp, 7);

    %------------------------------------------
    %% 6) Compute model-predicted values for PLOT
    %------------------------------------------
    if ~isempty(dateData_inc)
        midDaily_inc  = proportionToFitFunc(obj.fitProportions_inc, modelSolution , bothSim_inc);
        midSol_inc    = midDaily_inc(bothSim_inc);
        midSol_inc    = midSol_inc(:);
        t_inc         = t_inc(:);
        dateData_inc_common = dateData_inc_common(:);
    end

    if ~isempty(dateData_hosp)
        if obj.dataTypeLabel_hosp==2 && obj.hospital==1
            % Hospital prevalence from LOS-based estimate
            midDaily_hosp = midSol_hospPrev;
            midSol_hosp   = midDaily_hosp(bothSim_hosp);
        else
            % Daily hospitalizations or older model-based branch
            midDaily_hosp = proportionToFitFunc(obj.fitProportions_hosp, modelSolution, bothSim_hosp);
            midSol_hosp   = midDaily_hosp(bothSim_hosp);
        end

        midSol_hosp = midSol_hosp(:);
        t_hosp = t_hosp(:);
        dateData_hosp_common = dateData_hosp_common(:);
    end
    
    %------------------------------------------
    %% 7) Cut all dates and data at the Omicron wave start
    % We only fitted the extra days to more accurately calibrate that last few data points
    %------------------------------------------    
    if obj.cut_Omicron
        if ~isempty(dateData_inc)
            nInc = nnz(t_inc<obj.omicronWaveStart);

            midSol_inc=midSol_inc(1:nInc);
            sumData_inc=sumData_inc(1:nInc);
            smoothData_inc=smoothData_inc(1:nInc);
            dateData_inc_common=dateData_inc_common(1:nInc);
            t_inc=t_inc(1:nInc);
        end

        if ~isempty(dateData_hosp)
            nHosp = nnz(t_hosp<obj.omicronWaveStart);

            midSol_hosp=midSol_hosp(1:nHosp);
            sumData_hosp=sumData_hosp(1:nHosp);
            smoothData_hosp=smoothData_hosp(1:nHosp);  
            dateData_hosp_common=dateData_hosp_common(1:nHosp);
            t_hosp=t_hosp(1:nHosp);
        end
    end

    %------------------------------------------
    %% 7) Plot incidence
    %------------------------------------------
    figure
    
    if ~isempty(dateData_inc)
        subplot(4,4,[1 2 3 5 6 7])

        pl1 = plot(t_inc, midSol_inc, 'b:'); % Model-predicted incidence
        set(pl1, 'LineWidth', 1.5);
        hold on
    
        s1 = scatter(dateData_inc_common, sumData_inc, '+'); % Raw observed data points
        s1.MarkerEdgeColor = colourList(2);
    
        pa = plot(dateData_inc_common, smoothData_inc, '-'); % Smoothed daily data
        pa.Color = [colourList(2), .4];
    
        if obj.dataTypeLabel_inc==0
            ylabel("Daily Incidence")
        elseif obj.dataTypeLabel_inc==2
            ylabel("Incidence Prevalence")
        else
            ylabel("Daily Incidence")
        end

        vocLines(obj)
        legend("Model Prediction","Data","Smoothed Data",'Location','northwest')
        if obj.cut_Omicron
            XL = get(gca, 'XLim');
            xlim([XL(1) obj.omicronWaveStart])
        end
        ylim([0,5000])
    
        % Set 7 x-ticks for incidence plot
        XL = get(gca, 'XLim');
        dateRange = linspace(XL(1), XL(2), 7); % Creates 7 equally spaced ticks between the limits
        xticks(dateRange);
        datetick('x','mmm-yy','keepticks')
    end

    %------------------------------------------
    %% 8) Plot hospital
    %------------------------------------------
    if ~isempty(dateData_hosp)
        subplot(4,4,[9 10 11 13 14 15])
    
        pl2 = plot(t_hosp, midSol_hosp, 'b:'); % Model-predicted hospital
        set(pl2, 'LineWidth', 1.5);
        hold on
    
        s2 = scatter(dateData_hosp_common, sumData_hosp, '+'); % Raw observed data
        s2.MarkerEdgeColor = colourList(2);
    
        pb = plot(dateData_hosp_common, smoothData_hosp, '-'); % Smoothed daily data
        pb.Color = [colourList(2), .4];
    
        if obj.dataTypeLabel_hosp==0
            ylabel("Daily Hospitalizations")
        elseif obj.dataTypeLabel_hosp==2
            ylabel("Hospitalizations Prevalence")
        else
            ylabel("Daily Hospitalizations")
        end

        datetick('x','mmm-yy','keepticks')
        vocLines(obj)
        legend("Model Prediction","Data","Smoothed Data",'Location','northwest')    
        if obj.cut_Omicron
            XL = get(gca, 'XLim');
            xlim([XL(1) obj.omicronWaveStart])
        end
    
        % Set 7 x-ticks for hospital plot
        XL = get(gca, 'XLim');
        dateRange = linspace(XL(1), XL(2), 7); % Creates 7 equally spaced ticks between the limits
        xticks(dateRange);
        datetick('x','mmm-yy','keepticks')
    end

    %------------------------------------------
    %% 9) Title
    %------------------------------------------
    
    if obj.hospital == 0
        fitProb_text=strcat(", Inc=" ,num2str(obj.fitProportions_inc(2)),"$I_2+$",num2str(obj.fitProportions_inc(3)),"$I_3+$",num2str(obj.fitProportions_inc(4)),"$I_4$");
        sgtitle(strcat("Fitting cases only",fitProb_text),'Interpreter','latex')
    elseif obj.hospital == 1
        fitProb_text=strcat(", Hosp=" ,num2str(obj.fitProportions_hosp(4)),"$I_4$");
        sgtitle(strcat("Fitting hospitalizations only",fitProb_text),'Interpreter','latex')
    else
        fitProb_text=strcat(", Inc=" ,num2str(obj.fitProportions_inc(2)),"$I_2+$",num2str(obj.fitProportions_inc(3)),"$I_3+$",num2str(obj.fitProportions_inc(4)),"$I_4$",", Hosp=" ,num2str(obj.fitProportions_hosp(4)),"$I_4$");
        sgtitle(strcat("Fitting BOTH cases and hospitalizations",fitProb_text),'Interpreter','latex')
    end

    %% 10) Compute separate errors for each dataset
    if ~isempty(dateData_inc)&&~isempty(dateData_hosp)
        resid_len=min(length(smoothData_inc),length(midSol_inc));

        residInc=zeros(resid_len,1);
        for i=1:resid_len
            if smoothData_inc(i)~=0 && midSol_inc(i)~=0
                %residInc(i)=(smoothData_inc(i)-midSol_inc(i))./max(eps,midSol_inc(i));
                residInc(i)=(smoothData_inc(i)-midSol_inc(i))./max(eps,abs(smoothData_inc(i)+midSol_inc(i))/2); % percent error
            end
        end

        residIncSum=sum(abs(residInc));

        resid_len=min(length(smoothData_hosp),length(midSol_hosp));

        residHosp=zeros(resid_len,1);
        for i=1:resid_len
            if smoothData_hosp(i)~=0 && midSol_hosp(i)~=0
                %residHosp(i)=(smoothData_hosp(i)-midSol_hosp(i))./max(eps,midSol_hosp(i));
                residHosp(i)=abs(smoothData_hosp(i)-midSol_hosp(i))./max(eps,abs(smoothData_hosp(i)+midSol_hosp(i))/2);
            end
        end    
        residHospSum=sum(abs(residHosp));
        
        %% 12) Display final partial errors in a bar plot
        %------------------------------------------
        subplot(4,4,[8 12])
        bar(["Cases","Hospitalizations"], [residIncSum, residHospSum])
        ylabel("Sum of Percent Errors")
    end

    if obj.hospital==0
        residHospSum=NaN;
        resid_len=min(length(smoothData_inc),length(midSol_inc));
        residInc=zeros(resid_len,1);
        for i=1:resid_len
            if smoothData_inc(i)~=0 && midSol_inc(i)~=0
                residInc(i)=(smoothData_inc(i)-midSol_inc(i))./max(eps,abs(smoothData_inc(i)+midSol_inc(i))/2); % percent error
            end
        end
        residIncSum=sum(abs(residInc));

    elseif obj.hospital==1
        residIncSum=NaN;
        resid_len=min(length(smoothData_hosp),length(midSol_hosp));
        residHosp=zeros(resid_len,1);
        for i=1:resid_len
            if smoothData_hosp(i)~=0 && midSol_hosp(i)~=0
                residHosp(i)=abs(smoothData_hosp(i)-midSol_hosp(i))./max(eps,abs(smoothData_hosp(i)+midSol_hosp(i))/2);
            end
        end    
        residHospSum=sum(abs(residHosp));
    end

end


%% ----------------- Local Helper Functions ----------------- %% 

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
%   If dataIn has multiple columns, sum across them (assuming 2nd dimension is age).
%   Otherwise, return dataIn as is.
%
    if size(dataIn, 2) > 1
        mergedData = sum(dataIn, 2);
    else
        mergedData = dataIn;
    end
end