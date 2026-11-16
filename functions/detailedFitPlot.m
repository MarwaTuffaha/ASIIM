function detailedFitPlot(obj, pp, fitWindows, fast)
    arguments
        obj
        pp
        fitWindows = [0,0];
        fast = false;
    end

    % Process simulation output
    rp = postProcess(obj.yout);
    time = obj.time;

    % Ensure objective function value is finite
    % if ~isfinite(obj.objectiveFunctionValue)
    %     objValue = 1;
    % else
    %     objValue = obj.objectiveFunctionValue;
    % end

    % Extract case data depending on hospital or incidence dataset
    if obj.hospital==1
        dateData = obj.hospDates;
        data = obj.hospCases;
    else
        dateData = obj.incidenceDates;
        data = obj.incidenceCases;
    end

    % Validate case data and identify mitigation periods
    checkCaseData(data, obj);
    %mitList = pickMitagationWindows(obj,0,0);

    % Apply curvature penalty to smooth parameter estimates
    %ddk = obj.curvaturePenalty * diff(obj.kval,2).^2;
    %objValue = objValue + sum(ddk);

    % Convert date formats and align simulation time with observed data
    date0 = obj.dateList(1,:);
    dateList = datenum(obj.dateList(:,3), obj.dateList(:,2), obj.dateList(:,1));
    day0Simulation = datenum(date0(3), date0(2), date0(1));

    modelDates = time(1:end-1) + day0Simulation;
    
    if obj.hospital==1 && obj.dataTypeLabel_hosp==0
        modelDates = modelDates + obj.hospitalReportingLag;
    end
    
    % Align data and model using the actual model dates, not only the simulation date range
    bothData = ismember(dateData, modelDates);
    bothSim  = ismember(modelDates, dateData(bothData));
    
    t = modelDates(bothSim);

    %    (If hospital data is also daily or if you have different flags, adapt)
    % if obj.hospital == 1
    %     dataTypeLabel = obj.dataTypeLabel_hosp;  % hospital-only
    % else
    %     dataTypeLabel = obj.dataTypeLabel_inc;   % incidence or both
    % end

    if or(or(obj.hospital == 0 && obj.dataTypeLabel_inc==0,obj.hospital == 1 && obj.dataTypeLabel_hosp==0),and(obj.hospital == 2 && obj.dataTypeLabel_inc==0, obj.dataTypeLabel_hosp==0)) % data is daily
        % Convert to daily new cases/hospitalizations
        modelSolution = convertDaily(rp);
        % time is typically one less if we do daily differences
        time = time(1:end-1);

    elseif obj.hospital == 0 && obj.dataTypeLabel_inc==2  % We only have incidence prevalence data     
        modelSolution = zeros(size(rp(1).I,1), size(rp(1).I,2), 8);
        for i = 1:8
            modelSolution(:,:,i) = rp(i).I;
        end
        modelSolution = modelSolution(1:end,:,:);

    elseif obj.hospital == 1 && obj.dataTypeLabel_hosp==2  % We only have hospital prevalence data
        if obj.LOS_dist=="lognormal"
            [t_hosp_prev, midSol_hospPrev] = getPrevLognormal(obj, obj.LOS_Q1, obj.LOS_median, obj.LOS_Q3);
        elseif obj.LOS_dist=="PERT"
            [t_hosp_prev, midSol_hospPrev] = getPrevPERT(obj, obj.LOS_Q1, obj.LOS_median, obj.LOS_Q3);
        end
    else
        error("Code does not deal with these types/combinations of data. Either incidence and hospital daily, or only incidence prevalence, or only hospital prevalence")
    end

    % Compute model-predicted values
    if obj.hospital==1
        if obj.dataTypeLabel_hosp==0 %data is daily hosp
            midDaily = proportionToFitFunc(obj.fitProportions_hosp, modelSolution, bothSim);
            midSol = midDaily(bothSim);

        elseif obj.dataTypeLabel_hosp==2 %data is hosp prevalence

            % Important:
            % t_hosp_prev returned by getPrevLognormal is already absolute datenum.
            % Do NOT add day0Simulation.
            hospPrevDates = t_hosp_prev;

            dateData_round = round(dateData);
            hospPrevDates_round = round(hospPrevDates);

            % These masks are only for the fitted/model window.
            % They should NOT be used to plot all observed data points.
            bothData = ismember(dateData_round, hospPrevDates_round);
            bothSim  = ismember(hospPrevDates_round, dateData_round(bothData));

            t = hospPrevDates(bothSim);

            midDaily = midSol_hospPrev;
            midSol = midDaily(bothSim);

            % Force column vectors for plotting
            t = t(:);
            midSol = midSol(:);

        end        
    else
        midDaily = proportionToFitFunc(obj.fitProportions_inc, modelSolution, bothSim); %fitted model's prediction of daily incidence cases
        midSol = midDaily(bothSim); %Extracts values from midDaily only at time steps where observed data is available
    end

    % ---------------------------------------------------------------------
    % Data for the fitting/model window only
    % Used only for optional length checking.
    % ---------------------------------------------------------------------
    if size(data(bothData,:),2) > 1
        sumData = sum(data(bothData,:),2);
    else
        sumData = data(bothData,:);
    end    

    sumData = sumData(:);
    midSol = midSol(:);

    % Optional strict check for debugging
    if length(sumData) ~= length(midSol)
        error('Plot alignment problem: fitted data length = %d, model length = %d', ...
            length(sumData), length(midSol));
    end

    % ---------------------------------------------------------------------
    % Data for plotting all observed points
    % This is the key change: plot all data, not only data(bothData,:).
    % ---------------------------------------------------------------------
    if size(data,2) > 1
        sumDataFull = sum(data,2);
    else
        sumDataFull = data;
    end

    sumDataFull = sumDataFull(:);

    % Apply LOESS smoothing to all observed data
    smoothDataFull = smooth(sumDataFull, 7, 'loess');

    % Compute upper/lower bounds of estimated values
    %objScaled = (1/1) * objValue;
    %sumDaily = sum(modelSolution(:,:,3),2) + sum(modelSolution(:,:,4),2)+sum(modelSolution(:,:,7),2) + sum(modelSolution(:,:,8),2); % Upper bound I3+I4
    %daily = sum(modelSolution(:,:,4),2)+sum(modelSolution(:,:,8),2); % Lower bound I4
    %maxSol = sumDaily(bothSim); %Filters the upper bound to only include observed time steps 
    %minSol = daily(bothSim); %Filters the lower bound to only include observed time steps 

    %timeBreaks = timeBreaksFunc(obj.dateList);
    %lengthTimeBreaks = length(timeBreaks);
    %dateList = datenum(obj.dateList(:,3), obj.dateList(:,2), obj.dateList(:,1));
    %ddkOverObj = sum(ddk) / objValue;
    fontSize = 18;

    % Plot the estimated fit and observed data
    ax = gca;
    ax.FontSize = fontSize - 3;
    %plot(t, maxSol, 'c-.', 'LineWidth', 1.5); % Upper bound of estimates I3+I4
    hold on;

    pl = plot(t, midSol, 'b:'); % Model-predicted values in fitting/model window

    if ~isempty(pl)
        set(pl, 'LineWidth', 1.5);
    end

    %plot(t, minSol, 'k', 'LineWidth', 1.5); % Lower bound of estimates I3

    % Plot ALL observed data points
    s1 = scatter(dateData, sumDataFull, '+');
    s1.MarkerEdgeColor = colourList(2);

    % Plot smoothed ALL observed data
    pa = plot(dateData, smoothDataFull, '-');
    pa.Color = [colourList(2), .4];

    ylabel('Cases', 'interpreter', 'latex', 'FontSize', fontSize);
    xlabel('Time (Days)', 'interpreter', 'latex', 'FontSize', fontSize);
    %set(gca, 'YLim', [0, min(1.4 * max(sumData), 1.2 * max(maxSol))]);
    %set(gca,'YLim',[0,max(min(1.2*max(data),1.2*max(maxSol)),.1)]);
    set(gca,'YLim',[0,1.2*max(sumDataFull)]);
    
     if obj.hospital == 1
        verticalLines(obj, obj.hospDates);
    else
        verticalLines(obj, obj.incidenceDates);
     end
     obj.vocLines;

    if ~isequal(fitWindows, [0,0])
        xl(1) = xline(obj.dateNumList(fitWindows(1)), '-.');
        xl(2) = xline(obj.dateNumList(fitWindows(end)), '-.');
        xl(3) = xline(obj.dateNumList(fitWindows(1)), 'c-.');
        xl(4) = xline(obj.dateNumList(fitWindows(end)), 'c-.');
        for i = 1:4
            xl(i).LabelHorizontalAlignment = 'left';
            xl(i).LabelVerticalAlignment = 'middle';
            xl(i).FontSize = 8;
            xl(i).LineWidth = 1.5;
        end
        xl(1).Label = 'fit start';
        xl(2).Label = 'fit end';
    end     

     if obj.hospital == 1
        %legend('$I_3+I_4$', ['$' num2str(obj.fitProportions_hosp(4)) 'I_4$'], '$I_4$', 'hospitalization data', 'smoothed hospitalization data', 'interpreter', 'latex', 'FontSize', fontSize-3, 'Location', 'northwest');
        legend('Model Prediction', 'hospitalization data', 'smoothed hospitalization data', 'interpreter', 'latex', 'FontSize', fontSize-3, 'Location', 'northwest');
        if obj.dataTypeLabel_hosp==0 %daily data
            title(['Region: ', obj.region{1}, ', Simulation: ', obj.experimentName{1}, ', New Daily Cases']);
        elseif obj.dataTypeLabel_hosp==2 %prevalence data 
            title(['Region: ', obj.region{1}, ', Simulation: ', obj.experimentName{1}, ', Prevalence']);
        end
    else
        %legend('$I_3+I_4$', ['$' num2str(obj.fitProportions_inc(2)) 'I_2+' num2str(obj.fitProportions_inc(3)) 'I_3+' num2str(obj.fitProportions_inc(4)) 'I_4$'], '$I_4$','incidence data', 'smoothed incidence data', 'interpreter', 'latex', 'FontSize', fontSize-3, 'Location', 'northwest');
        legend('Model Prediction','incidence data', 'smoothed incidence data', 'interpreter', 'latex', 'FontSize', fontSize-3, 'Location', 'northwest');
        if obj.dataTypeLabel_inc==0 %daily data
            title(['Region: ', obj.region{1}, ', Simulation: ', obj.experimentName{1}, ', New Daily Cases']);
        elseif obj.dataTypeLabel_inc==2 %prevalence data 
            title(['Region: ', obj.region{1}, ', Simulation: ', obj.experimentName{1}, ', Prevalence']);
        end
     end

    datetick('x', 'mmm-yy', 'keepticks');
    
    h=gcf;
    set(h,'PaperType','uslegal');
    set(h,'PaperPositionMode','auto');
    set(h,'PaperOrientation','landscape');
    
    %set(h,'Position',[50 50 1200 800]);
    set(h,'Position',pp.position);    
end