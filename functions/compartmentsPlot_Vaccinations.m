function compartmentsPlot_Vaccinations(obj, BoosterLabel,Title,ageGroups)

arguments
    obj
    BoosterLabel=1; %To plot booster coverage (1), or not (0)
    Title="";
    ageGroups = 1:16;
end

% compartmentsPlot Plots compartment data as a percentage of the population
% and plots cumulative compartments
%
%
% Description:
%   This function generates two figures visualizing the population compartments 
%   (such as Susceptible, Infected, Exposed, Vaccinated) as a percentage of the total 
%   population for specified age groups. The compartments are plotted over time.
%
% Input Arguments:
%   obj       - Object containing the simulation data and configuration.
%   ageGroups - Vector specifying the age groups to include in the plots 
%               (optional, default is 1:16).
%
% Examples:
%   % Plot compartments for all age groups
%   [f1, f2] = compartmentsPlot(simulationObject);
%
%   % Plot compartments for specific age groups
%   [f1, f2] = compartmentsPlot(simulationObject, [1, 5, 10]);


r = obj.post;
tpop = sum(obj.population(ageGroups));
timeStart = obj.dateNumList(1);
time = obj.time + timeStart;

% Initialize compartment variables
len = length(sum(r(1).S(:, ageGroups), 2));
V1 = zeros(len, 1);
V2 = zeros(len, 1);
VB = zeros(len, 1);
CV1 = zeros(len, 1);
CV2 = zeros(len, 1);
CVB = zeros(len, 1);

% Aggregate compartments


for i = 1
    V1 = V1 + sum(r(i).V1(:, ageGroups), 2) + sum(r(i).V2(:, ageGroups), 2) + sum(r(i).V3(:, ageGroups), 2) + sum(r(i).V4(:, ageGroups), 2);
end
for i = 2
    V2 = V2 + sum(r(i).V1(:, ageGroups), 2) + sum(r(i).V2(:, ageGroups), 2) + sum(r(i).V3(:, ageGroups), 2) + sum(r(i).V4(:, ageGroups), 2);
end
VB = VB + sum(r(1).VB(:, ageGroups), 2) + sum(r(2).VB(:, ageGroups), 2) + sum(r(3).VB(:, ageGroups), 2) + sum(r(4).VB(:, ageGroups), 2);

for i = 1:4
    CV1 = CV1 + sum(r(i).CV1(:, ageGroups), 2);
    CV2 = CV2 + sum(r(i).CV2(:, ageGroups), 2);
    CVB = CVB + sum(r(i).CVB(:, ageGroups), 2);
end

font_size=15;

% Plot 1: Compartments as percent of population
f1 = figure('Name', 'Vaccination Percent of Population');
hold on
plot(time, V1/tpop, '-.','LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
plot(time, V2/tpop, '-.','LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980])
if BoosterLabel
    plot(time, VB/tpop, '-.','LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250])
end
if obj.hospital ==1
    verticalLines(obj, obj.hospDates, 'top')
else
    verticalLines(obj, obj.incidenceDates, 'top')
end
datetick('x', 'mmm-yy', 'keepticks')
legend( 'V1', 'V2', 'VB', 'Location', 'southwest')
savePlot(obj, 'vaccinated_compartments')
close

% Plot 2: Cumulative compartments
f2 = figure('Name', 'Cummulative Vaccinations');
hold on
if BoosterLabel
    vaccineCoveragePlotFunc(obj,1:3)
else
    vaccineCoveragePlotFunc(obj,1:2)
end   
plot(time, CV1/tpop,'LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
plot(time, CV2/tpop,'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980])
if BoosterLabel
    plot(time, CVB/tpop,'LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250])
end

line_width=2;
font_size=13;
xline(obj.alphaWaveStart,'-.','LineWidth',line_width,'Label',"Alpha","FontSize",font_size,'LabelHorizontalAlignment','left','Color',[0, 0.5, 0]);
xline(obj.deltaWaveStart,'-.','LineWidth',line_width,'Label',"Delta","FontSize",font_size,'Color',[0, 0.5, 0]);
xline(obj.vaccinationStart,'--','LineWidth',line_width,'Label',"Vacc. Start","FontSize",font_size,'LabelHorizontalAlignment','right')

if BoosterLabel
    xline(obj.vaccinationBoostStart,'--','LineWidth',line_width,'Label',"Booster Start","FontSize",font_size,'LabelHorizontalAlignment','right')
end
XL = get(gca, 'XLim');
xlim([738100 XL(2)]) % start plotting right before vaccination starts

if obj.cut_Omicron
    XL = get(gca, 'XLim');
    xlim([XL(1) obj.omicronWaveStart])
end

set(gca,'fontsize',font_size)

% Set 7 x-ticks for incidence plot
XL = get(gca, 'XLim');
dateRange = linspace(XL(1), XL(2), 8); % Creates 7 equally spaced ticks between the limits
xticks(dateRange);
datetick('x','mmm-yy','keepticks')

ylabel("Population Proportion",'Interpreter','latex')


if BoosterLabel
    leg=legend('Coverage Data - Dose1','Coverage Data - Dose2','Coverage Data - Dose3','Cumulative $V^1$', 'Cumulative $V^2$', 'Cumulative $V^B$', 'Location', 'southeast','Interpreter','latex');
else
    leg=legend('Coverage Data - Dose1','Coverage Data - Dose2','Cumulative $V^1$', 'Cumulative $V^2$', 'Location', 'southeast','Interpreter','latex');
end
leg.Position(3) = leg.Position(3) + 0.05; % widen box
savePlot(obj, 'Cumm_vacc')

end
