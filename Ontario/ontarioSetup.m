function sim = ontarioSetup()
%ONTARIOSIMULATIONDEFINITION Initialize a simulation for COVID-19 in Ontario.
% This function sets up a simulation environment using the asiim class,
% configuring parameters specific to COVID-19 dynamics in Ontario

%% Initialize an experiment object as the asiim class
sim = asiim; % Create an instance of the asiim class to represent the simulation.

%% Experiment settings
% Set paths and identifiers for the simulation data and results.
ontarioFolder = fileparts(mfilename('fullpath'));
sim.folder = ontarioFolder;
sim.dataFolder = fullfile(ontarioFolder, 'data');
sim.experimentName = 'ontario-1-preOm';
sim.description = 'COVID-19 FitHosp in Ontario, PreOmicron';
sim.shortTitle = 'ONT-FitHosp-PreOmic';

sim.cut_Omicron=1;

%% Create results folder
% Check if the results directory exists, and create it if it does not.
resultsPath = fullfile(sim.experimentDir,sim.folder);
if ~exist(resultsPath, 'dir')
    mkdir(resultsPath); % Create the directory if it does not exist.
    disp(['Created directory: ', resultsPath]);
end

% Calibration writes its fitting time to the results directory.
if ~exist(fullfile(sim.folder, 'results'), 'dir')
    mkdir(fullfile(sim.folder, 'results'));
end

%% Load date list
% Load the list of dates for the simulation from a text file.
sim.dateListHandle = fullfile(sim.dataFolder,'dateList.csv');
dateTable = readtable(sim.dateListHandle);
dateListDT = datetime(dateTable.Var1,'InputFormat','yyyy-MM-dd');
sim.dateList = [dateListDT.Day dateListDT.Month dateListDT.Year];

%% Set simulation stop date
% Define the calendar date when the simulation should stop.
%sim.calenderStopDate = datestr(739282); % This uses MATLAB's datenum format.
sim.calenderStopDate = datestr(sim.omicronWave); % This uses MATLAB's datenum format.
%% Load data

sim.hospital = 2; % Sets the flag to '2' to indicate using both incidence and hospitalization data for fitting
sim.incidenceFile = fullfile(sim.dataFolder,'cases_data.csv');
incidenceTable = readtable(sim.incidenceFile);
sim.incidenceDates=datenum(incidenceTable.date);
sim.incidenceCases=incidenceTable.value;
sim.dataTypeLabel_inc=1; %(1 means cummulative and 0 means daily)

if sim.dataTypeLabel_inc==1 % if cummulative, change to daily
    sim.incidenceCases=reshape(sim.incidenceCases,length(sim.incidenceCases),1);
    sim.incidenceCases = [sim.incidenceCases(1);diff(sim.incidenceCases)];  % Compute daily new cases
    sim.incidenceDates = sim.incidenceDates(1:end);  % Adjust time indexing
    sim.dataTypeLabel_inc=0;
end

% data = load(fullfile(sim.dataFolder,'hospital_data.mat'));
% sim.hospDates = data.data.dates;
% sim.hospCases = data.data.cases;

hospFile = fullfile(sim.dataFolder,'hospital_data_unified.csv');
hospTable = readtable(hospFile);
sim.hospDates=datenum(hospTable.date);
sim.hospCases=hospTable.hosp_daily;
sim.dataTypeLabel_hosp=0; %(1 means cummulative and 0 means daily)

if sim.dataTypeLabel_hosp==1 % if cummulative, change to daily
    sim.hospCases=reshape(sim.hospCases,length(sim.hospCases),1);
    sim.hospCases = [sim.hospCases(1);diff(sim.hospCases)];  % Compute daily new cases
    sim.hospDates = sim.hospDates(1:end);  % Adjust time indexing
    sim.dataTypeLabel_hosp=0;
end

%% Error check
% Check the case data for any inconsistencies or errors.
checkCaseData(sim.incidenceCases, sim);

%% Load population data
% Load population data from a file and store it in the simulation object.
sim.populationHandle = fullfile(sim.dataFolder ,'population.txt');
populationTable = readtable(sim.populationHandle);
sim.population = populationTable.Population';

%% Load k-value
% Load k-values, which are parameters used in the model, from a file.
sim.kValueHandle = fullfile(sim.dataFolder ,'kvalue.csv');
kvalTable = readtable(sim.kValueHandle);
sim.kval = kvalTable.Var1;
% if length(sim.kval) ~= length(sim.dateList)
%     disp("Error: kValue List not equal length of date list")
% end

%% Load demographics

sim.q = 0; % Set quarantine parameter, modify as needed.
sim.mu = zeros(1,16); % Set mortality rates, modify as needed.

%% Load omega data
% Load omega data, which are parameters affecting disease transmission.
sim.omegaHandle = fullfile(sim.dataFolder ,'omega.txt');
omegaTable = readtable(sim.omegaHandle);
omegaDateNum = datenum(datetime(omegaTable.date));
omegaIndex = ismember(omegaDateNum,sim.dateNumList);
fieldNames = {'omega', 'omega2', 'omega3', 'omega4'};
for i = 1:length(fieldNames)
    sim.(fieldNames{i}) = omegaTable.omega(omegaIndex);
end

%% Load pax data
% Load data on Paxlovid, a therapeutic agent used in the simulation.
sim.paxHandle = fullfile(sim.dataFolder, 'pax.txt');
paxTable = readtable(sim.paxHandle);
paxIndex = ismember(datenum(datetime(paxTable.date)), sim.dateNumList);
sim.pax = paxTable.pax(paxIndex);

%% Load contact modifiers
% Load data that modify contact rates in different settings.
contactModifersTable = readtable(fullfile(sim.dataFolder ,'contactModifiers.csv'));
sim.contactModifers = contactModifersTable{:,2:4};

%% Load contact matrices
% Load matrices defining contact rates in different settings like home, work, etc.
fileNames = {'contactHome.csv', 'contactWork.csv', 'contactOtherPlaces.csv', 'contactSchool.csv'};
fieldNames = {'H', 'W', 'O', 'S'};
for i = 1:length(fileNames)
    fullFileName = fullfile(sim.dataFolder, fileNames{i});
    sim.(fieldNames{i}) = readmatrix(fullFileName);
end

%% Load perturbation matrices
% Load matrices defining specific perturbations to contact matrices.
prefixes = {'SF', 'WF', 'OF'};
prefixNumbers = [1, 2, 3];
for i = 1:length(prefixes)
    for j = 1:length(prefixNumbers)
        fileName = ['contact', prefixes{i}, num2str(prefixNumbers(j)), '.txt'];
        fullFileName = fullfile(sim.dataFolder, fileName);
        sim.(sprintf('%s%d', prefixes{i}, prefixNumbers(j))) = load(fullFileName);
    end
end

%% Fit proportions
% Set the proportions for fitting the model to data.
sim.fitProportions_inc = [0 0 0.6 1]; %fitProportions for incidence
sim.fitProportionsLabel_inc = '0-0-0.6-1';

sim.fitProportions_hosp = [0 0 0 0.15]; %fitProportions for hospitalizations
sim.fitProportionsLabel_hosp = '0-0-0-.15';

%% Load vaccination data
% Load data related to vaccination rates and coverage.
sim.vaccinateHandle = fullfile(sim.dataFolder, 'vaccinate.csv');
vaccinateTable = readtable(sim.vaccinateHandle);
sim.vaccinate = vaccinateTable.vaccinate;

%% Load vaccination coverage
% Load detailed data on vaccination coverage over time.
sim.vaccineCoverageHandle = fullfile(sim.dataFolder ,'vaccineCoverage.csv');
sim.vaccineCoverageLabel=1;
sim.vacTargetPopulationHandle = fullfile(sim.dataFolder ,'vaccineTargetPopulation.txt');
vacTargetTable = readtable(sim.vacTargetPopulationHandle);
vacTargetData = table2array(vacTargetTable(:,2:end));
vacTargetDates = datenum(vacTargetTable.Var1);
sim.vacTargetPopulation = vacTargetData(ismember(vacTargetDates,sim.dateNumList),:);
sim.newVacTarget = vacTargetData(ismember(vacTargetDates,sim.dateNumList),:);

%%

% Cut all data after Omicron
% Keep some data points after Omicron to help fitting the last few points
sim.morePoints=50;
h=find(sim.hospDates==sim.omicronWave+sim.morePoints);
sim.hospDates=sim.hospDates(1:h);
sim.hospCases=sim.hospCases(1:h);

i=find(sim.incidenceDates==sim.omicronWave+sim.morePoints);
sim.incidenceDates=sim.incidenceDates(1:i);
sim.incidenceCases=sim.incidenceCases(1:i);

dateNumList = datenum(sim.dateList(:,[3,2,1]));
d=find(dateNumList>=sim.omicronWave+sim.morePoints,1,'first');
sim.dateList=sim.dateList(1:d,:);
sim.pax=sim.pax(1:d);
sim.kval=sim.kval(1:d);
sim.omega=sim.omega(1:d,:);

%sim.alphaInput = [1, 0.8, 0.4]; % Weaker protection
%sim.alphaInput = [1, 2/3, 1/3]; % Current Baseline
%sim.alphaInput = [1, 0.5, 0.25]; % Moderately stronger protection
%sim.alphaInput = [1, 1/3, 1/6]; % Stronger protection
%sim.alphaInput = [1, 0.2, 0.1]; % Much stronger protection


%% Model parameters

%% Fit properties
sim.curvaturePenaltyPercent = 1e-3;
sim.curvaturePenalty = 1e1;
sim.populationFactor = 6;
sim.timeLimit = 2;       
sim.windowLoopStart = 3;
sim.windowLoopEnd = 102;  
sim.windowSize = 16;
sim.windowStep = 1;     
sim.upperBound = 16; 
sim.lowerBound = 0.1;
sim.maxStallTime = 60*5; 
sim.maxGenerations = 240;

%% File handles
sim.forceFitHandle = fullfile(sim.dataFolder ,'forceFit.csv');
sim.exitFlagHandle = fullfile(sim.dataFolder ,'exitFlags.txt');

%% Set booster dose decay
sim.boosterDeacyFlag = false;

%% Run simulation
sim.simulate(); % Call the simulate method to run the simulation based on the configured parameters.

end
