function  compartmentsPlot_v3(obj, ageGroups)
% compartmentsPlot_v2 Plots detailed compartment data for specified scenarios and age groups
%
% Syntax:
%   compartmentsPlot_v3(obj)
%   compartmentsPlot_v3(obj, ageGroups)
%
% Description:
%   This function generates detailed plots of population compartments 
%   (such as Susceptible, Infected, Vaccinated, and Boosted) for multiple scenarios 
%   and age groups over time. It produces a tiled layout of subplots, each representing 
%   different compartments. The function allows for customization 
%   by specifying age groups and the number of scenarios to display.
%
% Input Arguments:
%   obj          - Object containing the simulation data and configuration.
%   ageGroups    - Vector specifying the age groups to include in the plots 
%                  (optional, default is 1:16).
%
% Examples:
%   % Plot compartments for all age groups
%   [f1, f2] = compartmentsPlot_v3(simulationObject);
%
%   % Plot compartments for specific age groups
%   [f1, f2] = compartmentsPlot_v3(simulationObject, [1, 5, 10]);


arguments
    obj
    ageGroups = 1:16;
end

r = obj.post;
tpop = sum(obj.population(ageGroups));
timeStart = obj.dateNumList(1);
time = obj.time + timeStart;

% Initialize compartment data
S = zeros(length(r(1).S), 4);
SV=zeros(length(r(1).S), 4);
I = zeros(length(r(1).I), 3);
IV = zeros(length(r(1).I), 3);
VB = zeros(length(r(1).VB), 4);
V1 = zeros(length(r(1).V1), 4);
V2 = zeros(length(r(1).V2), 4);

% Sum over age groups and scenarios
for i = 1:4
    S(:, i) = sum(r(i).S(:, ageGroups), 2);
    SV(:, i) = sum(r(i+4).S(:, ageGroups), 2);
    V1(:, i) = sum(r(1).(['V', num2str(i)])(:, ageGroups), 2);
    V2(:, i) = sum(r(2).(['V', num2str(i)])(:, ageGroups), 2);
    VB(:,i) = sum(r(i).VB(:,ageGroups),2);
end

for i = 2:4
    I(:, i-1) = sum(r(i).I(:, ageGroups), 2);
    IV(:, i-1) = sum(r(i+4).I(:, ageGroups), 2);
end

% Generate legend labels dynamically
legendLabelsS = arrayfun(@(n) sprintf('$S_{%d}$', n), 1:size(S,2), 'UniformOutput', false);
legendLabelsSV = arrayfun(@(n) sprintf('$S^V_{%d}$', n), 1:size(SV,2), 'UniformOutput', false);
legendLabelsI = arrayfun(@(n) sprintf('$I_{%d}$', n), (1:size(I,2))+1, 'UniformOutput', false);
legendLabelsIV = arrayfun(@(n) sprintf('$I^V_{%d}$', n), (1:size(IV,2))+1, 'UniformOutput', false);
legendLabelsV1 = arrayfun(@(n) sprintf('$V^1_{%d}$', n), 1:size(V1,2), 'UniformOutput', false);
legendLabelsV2 = arrayfun(@(n) sprintf('$V^2_{%d}$', n), 1:size(V2,2), 'UniformOutput', false);
legendLabelsVB = arrayfun(@(n) sprintf('$V^B_{%d}$', n), 1:size(VB,2), 'UniformOutput', false);

figure('Name', 'S compartments');
tiledlayout(1,2);

% S compartments
subPlotCompartment2(obj,time, S, tpop, legendLabelsS, '$S$ Compartments');
ylabel("Population Proportion",'Interpreter','latex')

% SV compartments
subPlotCompartment2(obj,time, SV, tpop, legendLabelsSV, '$S^V$ Compartments');

savePlot(obj, 'S compartments');


figure('Name', 'I compartments');
tiledlayout(1,2);

% I compartments
subPlotCompartment2(obj,time, I, tpop, legendLabelsI, '$I$ Compartments',1);
ylabel("Population Proportion",'Interpreter','latex')

% I compartments
subPlotCompartment2(obj,time, IV, tpop, legendLabelsIV, '$I^V$ Compartments',1);

savePlot(obj, 'I compartments');


figure('Name', 'V compartments');
tiledlayout(1,2);
% V1 compartments
subPlotCompartment2(obj,time, V1, tpop, legendLabelsV1, '$V^1$ Compartments');
ylabel("Population Proportion",'Interpreter','latex')

% V2 compartments
subPlotCompartment2(obj,time, V2, tpop, legendLabelsV2, '$V^2$ Compartments');

% % VB compartments
% subPlotCompartment2(obj,time, VB, tpop, legendLabelsVB, '$V^B$ Compartments');
% ylim([0,0.4])

savePlot(obj, 'V compartments');

end
