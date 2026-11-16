function vaccineCoveragePlotFunc(obj,doses,Coverage_or_Number)

arguments
    obj
    doses=1:3; % doses to plot coverage for
    Coverage_or_Number="c"; % Label for whether we plot the coverage as a percentage of the population ("c") or the absolute number of doses ("d")
end
%vaccineOverPlotFunc plot vaccine deltas
%   Detailed explanation goes here
%scale = max([max(obj.vaccineCoverage(:,1)),...
%    max(obj.vaccineCoverage(:,2)),...
%    max(obj.vaccineCoverage(:,1))]);
min_len=min(sum(obj.vaccinate),length(obj.vaccineCoverage(:,1)));
dates=obj.dateNumList(obj.vaccinate);
hold on

if ~isnumeric(doses)
    if doses=="second+booster"
        if Coverage_or_Number=="c"
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,1)/100, '-.','LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,2)/100+obj.vaccineCoverage(1:min_len,3)/100, '-.','LineWidth',1.5,'Color',[0.4940, 0.1840, 0.5560])
        elseif Coverage_or_Number=="d"
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,1)/100*obj.totalPopulation, '-.','LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,2)/100*obj.totalPopulation+obj.vaccineCoverage(1:min_len,3)/100*obj.totalPopulation, '-.','LineWidth',1.5,'Color',[0.4940, 0.1840, 0.5560])
        end
    end
else
    for d=doses
        if Coverage_or_Number=="c"
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,d)/100, '-.','LineWidth',1.5)
        elseif Coverage_or_Number=="d"
            plot(dates(1:min_len),obj.vaccineCoverage(1:min_len,d)/100*obj.totalPopulation, '-.','LineWidth',1.5)
        end

    end
end

end