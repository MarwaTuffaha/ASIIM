function subPlotCompartment2(obj,time, data, tpop, legendLabels, titleStr, InfectedFlag,xtickNumber)

arguments
    obj
    time
    data
    tpop
    legendLabels
    titleStr
    InfectedFlag=0;
    xtickNumber=6;
end
    nexttile;
    hold on;
    line_width=2;
    font_size=13;
    colors=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
    line_styles=["-","--",":","-."];
    
    for i = 1:size(data,2)
            plt=plot(time, data(:, i) / tpop,'LineWidth',line_width);
            if InfectedFlag==0
                plt.LineStyle=line_styles(i);
                plt.Color=colors(i,:);
            else
                plt.LineStyle=line_styles(i+1);
                plt.Color=colors(i+1,:);
            end
    end
   
  
    xline(obj.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
    xline(obj.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'Color',[0, 0.5, 0]);
    xline(obj.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')
    
    
    if ~obj.cut_Omicron
        xline(obj.omicronWaveStart,'-.','LineWidth',line_width,'Label',"Omicron","FontSize",font_size,'LabelHorizontalAlignment','right','Color',[0, 0.5, 0]);
        xline(obj.vaccinationBoostStart,'--','LineWidth',line_width,'Label',"Booster Start","FontSize",font_size,'LabelHorizontalAlignment','right')
    end        

    legend(legendLabels, 'Location', 'northwest', 'Interpreter', 'latex');
    
    if obj.cut_Omicron==1
        xlim([time(1), obj.omicronWaveStart])
    else
        xlim([time(1), time(end)])
    end

    % Convert dateNumList to datetime
    dateTimes = datetime(obj.dateNumList, 'ConvertFrom', 'datenum');

    for i=1:length(dateTimes)
        if obj.cut_Omicron==1 && obj.dateNumList(i)>=obj.omicronWaveStart
            dateTimes=dateTimes(1:i-1);
            break
        end
    end

    % Select approximately 3-4 evenly spaced dates
    numTicks = xtickNumber; % Adjust this number to control how many dates are shown
    tickIndices = round(linspace(1, length(dateTimes), numTicks)); % Evenly spaced indices
    
    % Ensure first tick is included
    tickIndices(1) = 1; % Force the first tick to be at the start
    selectedTicks = obj.dateNumList(tickIndices); % Corresponding datenums
    
    % Set X-axis ticks and labels
    xticks(selectedTicks);
    xticklabels(datestr(datetime(selectedTicks, 'ConvertFrom', 'datenum'), 'mmm-yy'));

    %datetick('x', 'mmm-yy', 'keepticks');
    title(titleStr,'Interpreter','latex');
    ax=gca;
    ax.FontSize = font_size;
end