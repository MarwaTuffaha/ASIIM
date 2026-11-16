function [t_inc, inc, t_hosp, hosp] = getPrediction(obj)

    rp = postProcess(obj.yout);
    %time = [obj.time(1):obj.dateNumList(end)-obj.dateNumList(1)]';
    time = obj.time;

    % Convert cummulative compartment to daily new cases/hospitalizations
    modelSolution = convertDaily(rp);
    % time is typically one less if we do daily differences
    time = time(1:end-1);


    date0 = obj.dateList(1,:);
    day0Simulation_inc = datenum(date0(3), date0(2), date0(1));
    day0Simulation_hosp = datenum(date0(3), date0(2), date0(1));

    t_inc = time + day0Simulation_inc;

    t_hosp = time + day0Simulation_hosp + obj.hospitalReportingLag;    



    % % Compute model-predicted values
    % if obj.hospital==0
    %     inc = proportionToFitFunc(obj.fitProportions_inc, modelSolution); %fitted model’s prediction of daily incidence cases
    %     inc = sum(inc,2);
    %     hosp=[]; 
    % end
    % 
    % if obj.hospital==1
    %     hosp = proportionToFitFunc(obj.fitProportions_hosp, modelSolution); %fitted model’s prediction of daily hospitalization cases
    %     hosp = sum(hosp,2);
    %     inc=[];
    % end
    % 
    % if obj.hospital==2


    % Compute model-predicted incidence only if incidence fit proportions exist
    if isempty(obj.fitProportions_inc)
        inc = [];
    else
        inc = proportionToFitFunc(obj.fitProportions_inc, modelSolution);
        inc = sum(inc,2);
    end
    
    % Compute model-predicted hospital output only if hospital fit proportions exist
    if isempty(obj.fitProportions_hosp)
        hosp = [];
    else
        hosp = proportionToFitFunc(obj.fitProportions_hosp, modelSolution);
        hosp = sum(hosp,2);
    end    
    %     inc = proportionToFitFunc(obj.fitProportions_inc, modelSolution); %fitted model's prediction of daily incidence cases
    %     inc = sum(inc,2);
    %     hosp = proportionToFitFunc(obj.fitProportions_hosp, modelSolution); %fitted model's prediction of daily hospitalization cases
    %     hosp = sum(hosp,2);
    % % end    

    if obj.cut_Omicron==1
        om_ind_inc=find(t_inc==obj.omicronWave,1,"first");
        if ~isempty(inc)
            t_inc=t_inc(1:om_ind_inc);
            inc=inc(1:om_ind_inc);
        end
        om_ind_hosp=find(t_hosp==obj.omicronWave,1,"first");
        if ~isempty(hosp)
            t_hosp=t_hosp(1:om_ind_hosp);
            hosp=hosp(1:om_ind_hosp);
        end
    end

    
end