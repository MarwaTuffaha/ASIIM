function [] = resutsPlotsSeroprevalenceImmunityBoost(sSub,rSub,data,plotProp,Title)
arguments
    sSub
    rSub
    data
    plotProp
    Title=sSub.experimentName;
end

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fileName = plotProp.fileName;
ageRange = plotProp.ageRange;
visible = plotProp.visible;
showTitle = plotProp.showTitle;
save = plotProp.save;
subPlot = plotProp.subPlot;
showLegend = plotProp.showLegend;

dateDataMax = max(data.dates(end));

if rem(length(ageRange),2)~=0
    disp('Age Range Error!')
end

simIndex = 1;

dateList = datetime(sSub(1).dateList(:,[3,2,1]));

time = rSub(simIndex).time+datenum(dateList(1));
rp = postProcess(rSub(simIndex).yout);
prev = inf(length(rp(1).S(:,1)),length(ageRange)/2,3);
ageIndex = 0;

% Write to simulation Log file

logFile = fopen(sSub(1).simulationLogFilePath, 'a');
if logFile == -1
    error('Unable to open the log file.');
end
% Write the debug information to the file
fprintf(logFile, '%s - resultsPlotSeroprevalenceBoost.m\n', datetime('now'));
fprintf(logFile, '- Length of time: %d\n', length(time));
fprintf(logFile, '- Last simulation time: %s (check data\\dateList)\n', datestr(time(end)));
fprintf(logFile, '- Last date in data: %s\n', datestr(data.dates(end)));

% Close the log file
fclose(logFile);

for j = ageRange(1):2:ageRange(end)-1
    ageIndex = ageIndex + 1;
    prev(:,ageIndex,3) = sum(rp(2).S(:,j:j+1),2) + ...
                            sum(rp(1).V1(:,j:j+1),2)+...
                            sum(rp(6).S(:,j:j+1),2);
    prev(:,ageIndex,2) = sum(rp(3).S(:,j:j+1),2) + ...
                            sum(rp(1).V2(:,j:j+1),2) + ...
                            sum(rp(7).S(:,j:j+1),2);
    prev(:,ageIndex,1) = sum(rp(4).S(:,j:j+1),2) + sum(rp(2).V1(:,j:j+1),2)...
                            +sum(rp(2).V2(:,j:j+1),2)+sum(rp(1).V3(:,j:j+1),2)...
                            +sum(rp(2).V3(:,j:j+1),2)+sum(rp(1).V4(:,j:j+1),2)...
                            +sum(rp(2).V4(:,j:j+1),2)...
                            +sum(rp(1).VB(:,j:j+1),2)+sum(rp(2).VB(:,j:j+1),2)...
                            +sum(rp(3).VB(:,j:j+1),2)+sum(rp(4).VB(:,j:j+1),2)...
                            +sum(rp(8).S(:,j:j+1),2);
    prev(:,ageIndex,4) = sum(rp(1).S(:,j:j+1),2)+sum(rp(5).S(:,j:j+1),2);
end

%     prev(:,ageIndex,4) = sum(rp(6).S(:,j:j+1),2)+...
%         sum(rp(7).S(:,j:j+1),2)+sum(rp(8).S(:,j:j+1),2);

% ageIndex = 0;
% for j = ageRange(1):2:ageRange(end)-1
% ageIndex = ageIndex + 1;
% prev(:,ageIndex,1) = sum(rp(1).V1(:,j:j+1),2)+sum(rp(1).V2(:,j:j+1),2)+...
%         sum(rp(1).V3(:,j:j+1),2)+sum(rp(1).V4(:,j:j+1),2);
% prev(:,ageIndex,2) = sum(rp(2).V1(:,j:j+1),2)+sum(rp(2).V2(:,j:j+1),2)+...
%         sum(rp(2).V3(:,j:j+1),2)+sum(rp(2).V4(:,j:j+1),2);
% end
    
if plotProp.percent
    ageIndex = 0;
    for j = ageRange(1):2:ageRange(end)-1
    ageIndex = ageIndex + 1;
    comPop(ageIndex) = sum(sSub(simIndex).population(j:j+1),2);
    end
    for i = 1:4
    prev(:,:,i) = prev(:,:,i)./(sum(sSub(simIndex).population(ageRange),2));
    end
    
    legendText = {...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 0-9',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 10-19',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 20-29',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 30-39',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 40-49',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 50-59',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 60-69',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2+ \sum_{i=1}^{4} V_i^B$, Age 70+',...
        '$S_3 + S_3^V + V_2^1$, Age 0-9',...
        '$S_3 + S_3^V + V_2^1$, Age 10-19',...
        '$S_3 + S_3^V + V_2^1$, Age 20-29',...
        '$S_3 + S_3^V + V_2^1$, Age 30-39',...
        '$S_3 + S_3^V + V_2^1$, Age 40-49',...
        '$S_3 + S_3^V + V_2^1$, Age 50-59',...
        '$S_3 + S_3^V + V_2^1$, Age 60-69',...
        '$S_3 + S_3^V + V_2^1$, Age 70+',...
        '$S_2 + S_2^V + V_1^1$, Age 0-9',...
        '$S_2 + S_2^V + V_1^1$, Age 10-19',...
        '$S_2 + S_2^V + V_1^1$, Age 20-29',...
        '$S_2 + S_2^V + V_1^1$, Age 30-39',...
        '$S_2 + S_2^V + V_1^1$, Age 40-49',...
        '$S_2 + S_2^V + V_1^1$, Age 50-59',...
        '$S_2 + S_2^V + V_1^1$, Age 60-69',...
        '$S_2 + S_2^V + V_1^1$, Age 70+'...
        '$S_1 + S_1^V$, Age 0-9',...
        '$S_1 + S_1^V$, Age 10-19',...
        '$S_1 + S_1^V$, Age 20-29',...
        '$S_1 + S_1^V$, Age 30-39',...
        '$S_1 + S_1^V$, Age 40-49',...
        '$S_1 + S_1^V$, Age 50-59',...
        '$S_1 + S_1^V$, Age 60-69',...
        '$S_1 + S_1^V$, Age 70+',...
        %'vaccination start',...
        %'booster vaccination start',...
        %'past/prediction'
        };
else
    legendText = {...
              '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 0-9',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 10-19',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 20-29',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 30-39',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 40-49',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 50-59',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 60-69',...
        '$S_4 + S_4^V + \sum_{i=3}^{4} V_i^1 + \sum_{i=1}^{4} V_i^2$, Age 70+',...
        '$S_3 + S_3^V + V_2^1$, Age 0-9',...
        '$S_3 + S_3^V + V_2^1$, Age 10-19',...
        '$S_3 + S_3^V + V_2^1$, Age 20-29',...
        '$S_3 + S_3^V + V_2^1$, Age 30-39',...
        '$S_3 + S_3^V + V_2^1$, Age 40-49',...
        '$S_3 + S_3^V + V_2^1$, Age 50-59',...
        '$S_3 + S_3^V + V_2^1$, Age 60-69',...
        '$S_3 + S_3^V + V_2^1$, Age 70+',...
        '$S_2 + S_2^V + V_1^1$, Age 0-9',...
        '$S_2 + S_2^V + V_1^1$, Age 10-19',...
        '$S_2 + S_2^V + V_1^1$, Age 20-29',...
        '$S_2 + S_2^V + V_1^1$, Age 30-39',...
        '$$S_2 + S_2^V + V_1^1$, Age 40-49',...
        '$$S_2 + S_2^V + V_1^1$, Age 50-59',...
        '$S_2 + S_2^V + V_1^1$, Age 60-69',...
        '$S_2 + S_2^V + V_1^1$, Age 70+',...
        '$S_1 + S_1^V$, Age 0-9',...
        '$S_1 + S_1^V$, Age 10-19',...
        '$S_1 + S_1^V$, Age 20-29',...
        '$S_1 + S_1^V$, Age 30-39',...
        '$S_1 + S_1^V$, Age 40-49',...
        '$S_1 + S_1^V$, Age 50-59',...
        '$S_1 + S_1^V$, Age 60-69',...
        '$S_1 + S_1^V$, Age 70+',...
        'Total Population'};
end

fontSize = plotProp.fontSize;
lineWidth = 1.5;
if ~subPlot
figure('Visible', visible)
end
        h=gcf;
        set(h,'PaperType','uslegal');
        set(h,'PaperPositionMode','auto');
        set(h,'PaperOrientation','landscape');
        set(h,'Position',[1 1 1800/1.5 1000/1.5]);
hold on
box on
%pl(1) = bar(time,sum(prev(:,:,2:4),3),'stacked');
b = bar(time,reshape(prev(:,:,[1 2 3 4]),size(prev(:,:,1),1),[],1),'stacked',...
    'BarWidth', 1,'FaceColor','flat');
len = size(prev(:,:,1),2);
alpha = linspace(.1,1,len);
for k = 1:len
    b(k).CData = [0.4660, 0.6740, 0.1880];
    %b(k).CData = [0, 0.4470, 0.741];
    b(k).FaceAlpha = alpha(k);
end
for k = len+1:2*len
    b(k).CData = 	[0.9290, 0.6940, 0.1250];
    %b(k).CData = 	[0.5, 0.3, 0.5];
    b(k).FaceAlpha = alpha(k-len);
end
for k = 2*len+1:3*len
    b(k).CData = [0.3010, 0.7450, 0.9330];
    %b(k).CData = [0.1, 0, 0.5];
    b(k).FaceAlpha = alpha(k-2*len);
end
for k = 3*len+1:4*len
    b(k).CData = [1, .2, .2];
    %b(k).CData = [0, 0, 0];
    b(k).FaceAlpha = alpha(k-3*len);
end

%     p(5) = xline(dateDataMax);
%     p(5).Color = [0,0,0,.5];
%     p(5).LineWidth = lineWidth-.5;
%     p(5).LineStyle = '-.';

line_width=2;
font_size=15;



if plotProp.verticalLinesFlag
    xline(sSub.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
    xline(sSub.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'Color',[0, 0.5, 0]);
    xline(sSub.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')
end

if showLegend
legend(b(end:-1:1),legendText(end:-1:1),...
    'interpreter','latex','location','northeastoutside','fontsize',fontSize-7);
end

set(gca,'FontSize',fontSize-4)
set(gca,'Ylim',[0,1])


if sSub.cut_Omicron==1
    xlim([time(1), sSub.omicronWaveStart])
end

% Convert dateNumList to datetime
dateTimes = datetime(sSub.dateNumList, 'ConvertFrom', 'datenum');

for i=1:length(dateTimes)
    if sSub.cut_Omicron==1 && sSub.dateNumList(i)>=sSub.omicronWaveStart
        dateTimes=dateTimes(1:i-1);
        break
    end
end

% Select approximately 3-4 evenly spaced dates
numTicks = 8; % Adjust this number to control how many dates are shown
tickIndices = round(linspace(1, length(dateTimes), numTicks)); % Evenly spaced indices

% Ensure first tick is included
tickIndices(1) = 1; % Force the first tick to be at the start
selectedTicks = sSub.dateNumList(tickIndices); % Corresponding datenums

% Set X-axis ticks and labels
xticks(selectedTicks);
xticklabels(datestr(datetime(selectedTicks, 'ConvertFrom', 'datenum'), 'mmm-yy'));

% title(['Sim-',num2str(sSub(simIndex).simNum)]);

if ~subPlot
%xlabel('time')
if ~plotProp.percent
    ylabel('population')
else
    ylabel('% population')
end

%title(strcat("Immunity Distribution, ",Title))
if showTitle
    if strcmp(sSub(simIndex).omegaLabel{1},'0')
        title(['Sim-',num2str(sSub(simIndex).simNum),', Vaccine type ',num2str(sSub(simIndex).vaccineTypeLabel),...
            ', no waning'])
    else
        title(['Sim-',num2str(sSub(simIndex).simNum),', Vaccine type ',num2str(sSub(simIndex).vaccineTypeLabel),...
            ', waning over ',sSub(simIndex).omegaLabel{1},' year per compartment'])
    end
end
end

if save
            %Create Results Folder
            if ~exist([sSub(1).folder,'/results/',date], 'dir')
                mkdir([sSub(1).folder,'/results/',date])
            end
            fullFileName = fullfile(sSub(1).folder,'/results/',date, fileName);
            %Save Figure
            print([getExperimentDir(),'\',fullFileName],plotProp.fileType,'-r300')
end
            
end

