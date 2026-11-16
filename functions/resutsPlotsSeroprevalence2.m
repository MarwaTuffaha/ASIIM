function [] = resutsPlotsSeroprevalence2(obj,ageRange)
% Here we plot seroprevalence to compare with data from Ontario


font_size=15;
line_width=2;

simIndex = 1;

ageRange_all=1:obj.N;

dateList = datetime(obj(1).dateList(:,[3,2,1]));

time = obj(simIndex).time+datenum(dateList(1));
rp = postProcess(obj(simIndex).yout);
antiN1 = zeros(length(rp(1).S(:,1)),length(ageRange_all)/2); % Will include S^V here
antiN1_prop = zeros(length(rp(1).S(:,1)),length(ageRange_all)/2);
antiN2 = zeros(length(rp(1).S(:,1)),length(ageRange_all)/2); % Will NOT include S^V here
antiN2_prop = zeros(length(rp(1).S(:,1)),length(ageRange_all)/2);
antiS = zeros(length(rp(1).S(:,1)),length(ageRange_all)/2);
antiS_prop=zeros(length(rp(1).S(:,1)),length(ageRange_all)/2);
ageIndex = 0;
for j = ageRange_all(1):2:ageRange_all(end)-1
    ageIndex = ageIndex + 1;
    antiN1(:,ageIndex) = sum(rp(2).S(:,j:j+1),2)+...
        sum(rp(3).S(:,j:j+1),2)+sum(rp(4).S(:,j:j+1),2)+...
        sum(rp(6).S(:,j:j+1),2)+...
        sum(rp(7).S(:,j:j+1),2)+sum(rp(8).S(:,j:j+1),2);
    antiN1_prop(:,ageIndex)=antiN1(:,ageIndex)/sum(obj.population(j:j+1));

    antiN2(:,ageIndex) = sum(rp(2).S(:,j:j+1),2)+...
        sum(rp(3).S(:,j:j+1),2)+sum(rp(4).S(:,j:j+1),2);
    antiN2_prop(:,ageIndex)=antiN2(:,ageIndex)/sum(obj.population(j:j+1));

    antiS(:,ageIndex)=sum(rp(2).S(:,j:j+1),2)+...
        sum(rp(3).S(:,j:j+1),2)+sum(rp(4).S(:,j:j+1),2)+...
        sum(rp(5).S(:,j:j+1),2)+sum(rp(6).S(:,j:j+1),2)+...
        sum(rp(7).S(:,j:j+1),2)+sum(rp(8).S(:,j:j+1),2)+...
        sum(rp(1).V1(:,j:j+1),2)+sum(rp(1).V2(:,j:j+1),2)+...
         sum(rp(1).V3(:,j:j+1),2)+sum(rp(1).V4(:,j:j+1),2)+...
         sum(rp(2).V1(:,j:j+1),2)+sum(rp(2).V2(:,j:j+1),2)+...
        sum(rp(2).V3(:,j:j+1),2)+sum(rp(2).V4(:,j:j+1),2)+...
        sum(rp(1).VB(:,j:j+1),2)+sum(rp(2).VB(:,j:j+1),2)+...
        sum(rp(3).VB(:,j:j+1),2)+sum(rp(4).VB(:,j:j+1),2);
    antiS_prop(:,ageIndex)=antiS(:,ageIndex)/sum(obj.population(j:j+1));


end
antiN1_total=sum(rp(2).S(:,ageRange),2)+...
        sum(rp(3).S(:,ageRange),2)+sum(rp(4).S(:,ageRange),2)+...
        sum(rp(6).S(:,ageRange),2)+...
        sum(rp(7).S(:,ageRange),2)+sum(rp(8).S(:,ageRange),2);
antiN1_prop_total=antiN1_total/sum(obj.population(ageRange));

antiN2_total=sum(rp(2).S(:,ageRange),2)+...
        sum(rp(3).S(:,ageRange),2)+sum(rp(4).S(:,ageRange),2);
antiN2_prop_total=antiN2_total/sum(obj.population(ageRange));

antiS_total=sum(rp(2).S(:,ageRange),2)+...
        sum(rp(3).S(:,ageRange),2)+sum(rp(4).S(:,ageRange),2)+...
        sum(rp(5).S(:,ageRange),2)+sum(rp(6).S(:,ageRange),2)+...
        sum(rp(7).S(:,ageRange),2)+sum(rp(8).S(:,ageRange),2)+...
        sum(rp(1).V1(:,ageRange),2)+sum(rp(1).V2(:,ageRange),2)+...
         sum(rp(1).V3(:,ageRange),2)+sum(rp(1).V4(:,ageRange),2)+...
         sum(rp(2).V1(:,ageRange),2)+sum(rp(2).V2(:,ageRange),2)+...
        sum(rp(2).V3(:,ageRange),2)+sum(rp(2).V4(:,ageRange),2)+...
        sum(rp(1).VB(:,ageRange),2)+sum(rp(2).VB(:,ageRange),2)+...
        sum(rp(3).VB(:,ageRange),2)+sum(rp(4).VB(:,ageRange),2);
antiS_prop_total=antiS_total/sum(obj.population(ageRange));


% import sero data (Ontario here)
addpath('C:\Users\Marwa T\Covid-ASIIM\plusCovid-ASIIM\sero_data_Ontario')
[Wmean,Wlow,Whigh] =sero_data_extraction();


% generate a single-color gradient from light to dark blue for N ages
N=obj.N;
%blueScale = linspace(1, 0.3, N/2)';
%cmap = [zeros(N/2,1), zeros(N/2,1), blueScale];

cmap = [
    0.1216 0.4667 0.7059   % blue
    1.0000 0.4980 0.0549   % orange
    0.1725 0.6275 0.1725   % green
    0.8392 0.1529 0.1569   % red
    0.5804 0.4039 0.7412   % purple
    0.5490 0.3373 0.2941   % brown
    0.8902 0.4667 0.7608   % pink
    0.4980 0.4980 0.4980   % gray
];





figure

p=[];
for i=1:N/2
    p(i)=plot(time,antiN1_prop(:,i), 'Color', cmap(i,:),'LineWidth',line_width);
    hold on
end


xline(obj.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')
    
title("$S_2+S_3+S_4+CIV_2+CIV_3+CIV_4$ Normalized per age",'Interpreter','latex')
ylabel("Proportion")

Legends=["Age 0-9","Age 10-19","Age 20-29","Age 30-39","Age 40-49","Age 50-59","Age 60-69","Age 70+"];
legend(p,Legends,'Location','northwest')

if obj.cut_Omicron
    XL = get(gca, 'XLim');
    xlim([XL(1) obj.omicronWaveStart])
end

% Set 7 x-ticks for incidence plot
XL = get(gca, 'XLim');
dateRange = linspace(XL(1), XL(2), 7); % Creates 7 equally spaced ticks between the limits
xticks(dateRange);
datetick('x','mmm-yy','keepticks')

set(gca,'fontsize',font_size)


figure
p=[];
for i=1:N/2
    p(i)=plot(time,antiS_prop(:,i), 'Color', cmap(i,:),'LineWidth',line_width);
    hold on
end



xline(obj.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')
    
title("$\sum_2^4(S)+\sum_1^4(S^V)+V$ Normalized per age",'Interpreter','latex')
ylabel("Proportion")

Legends=["Age 0-9","Age 10-19","Age 20-29","Age 30-39","Age 40-49","Age 50-59","Age 60-69","Age 70+"];
legend(p,Legends,'Location','northwest')

if obj.cut_Omicron
    xlim([time(1) obj.omicronWaveStart])
end

% Set 7 x-ticks for incidence plot
XL = get(gca, 'XLim');
dateRange = linspace(XL(1), XL(2), 7); % Creates 7 equally spaced ticks between the limits
xticks(dateRange);
datetick('x','mmm-yy','keepticks')

set(gca,'fontsize',font_size)

cmap=[0, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
    0.6350, 0.0780, 0.1840];

figure

hold on
%pl(1) = bar(time,sum(prev(:,:,2:4),3),'stacked');

pl(1)=plot(time,antiS_prop_total,'LineWidth',line_width+1,'Color',cmap(1,:));
pl(3)=plot(time,antiN1_prop_total,'LineWidth',line_width+0.5,'Color',cmap(2,:));
pl(4)=plot(time,antiN2_prop_total,'LineWidth',line_width,'Color',cmap(3,:));

xline(obj.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')
    
FirstAges=["0","5","10","15","20","25","30","35","40","45","50","55","60","65","70","75"];
title(strcat("Total, Ages=",num2str(FirstAges(ageRange(1))),"+"),'Interpreter','latex')
ylabel("Population Proportion")

if obj.cut_Omicron
    xlim([datenum(dateList(1)) obj.omicronWaveStart])
else
    XL = get(gca, 'XLim');
    xlim([datenum(dateList(1)) XL(2)])
end


XL = get(gca, 'XLim');

series = setdiff(Wmean.Properties.VariableNames,'week_end','stable');  % e.g., 'Anti-N estimate','Anti-S estimate'
%co = get(gca,'ColorOrder');

lw_adjust=[0,0.5]; %line_width adjustment


cmap=[0.8500, 0.3250, 0.0980;
    0, 0.4470, 0.7410];


for k = 1:numel(series)
    x   = datenum(Wmean.week_end);
    y   = Wmean.(series{k});    % already in percent
    ylo = Wlow.(series{k});
    yhi = Whigh.(series{k});

    c = cmap(k,:);

    % 95% CI band behind the line
    fill([x; flipud(x)], [ylo; flipud(yhi)], c, ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none');

    %Mean line
    if k==1
        pl_ind=5;
    elseif k==2
        pl_ind=2;
    end
    pl(pl_ind)=plot(x, y, 'LineWidth', 1.6+lw_adjust(k), 'Color', c,'LineStyle','--');
end

% Set 7 x-ticks for incidence plot
dateRange = linspace(XL(1), XL(2), 8); % Creates 7 equally spaced ticks between the limits
xticks(dateRange);
datetick('x','mmm-yy','keepticks')

%Legends=["$\sum_2^4(S)+\sum_1^4(S^V)+V$",'Anti-S data',"$\sum_2^4(S)+\sum_2^4S^V$","$\sum_2^4(S)$",'Anti-N data'];
Legends=["Anti-S model","Anti-S data",'Anti-N model upper bound','Anti-N model lower bound','Anti-N data'];


leg=legend(pl,Legends,'Interpreter','latex','Location','northwest');
%leg.ItemTokenSize = [50,50];   % gives more spacing
leg.Position(3) = leg.Position(3) + 0.1; % widen box

set(gca,'fontsize',font_size)
            
end

