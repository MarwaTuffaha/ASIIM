function sim = ItalySetup()

%% Initialize an experiment object as the assim class
sim = asiim; 

%% Experiment settings
sim.folder = fileparts(mfilename('fullpath'));
sim.dataFolder = fullfile(sim.folder, 'data');
if ~exist(fullfile(sim.folder, 'results'), 'dir')
    mkdir(fullfile(sim.folder, 'results'));
end
sim.experimentName = 'Italy';
sim.description = 'Italy';
sim.shortTitle = 'Italy';

%% Load date list
sim.dateListHandle = fullfile(sim.dataFolder, 'dateList.txt');

dateTable = readtable(sim.dateListHandle);
dateListDT = datetime(dateTable.Var1,'InputFormat','yyyy-MM-dd');
sim.dateList = [dateListDT.Day dateListDT.Month dateListDT.Year];

%% Set simulation stop date
%sim.calenderStopDate = datestr(739035);%datestr(738502+365); %until omicron (December 12 2021)
sim.calenderStopDate = datestr(sim.dateNumList(end));%datestr(738502+365); %until omicron (December 12 2021)
sim.cut_Omicron=1;
%% Load incidence data

rawData = readtable(fullfile(sim.dataFolder, 'incidenceData.txt'));
incDates = datetime(rawData.incidenceDates, 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Convert date strings to datetime objects
sim.incidenceDates = datenum(incDates);
sim.incidenceCases = rawData.newDailyCases;

%% Load hospital admissions data

rawData = readtable(fullfile(sim.dataFolder, 'OWD_hospAdmissionData.csv'));
hospDates = datetime(rawData.date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Convert date strings to datetime objects
sim.hospDates = datenum(hospDates);
sim.hospCases = rawData.smoothed_hospAdmissions;
sim.hospCasesSmoothLabel=1; % Data is already smoothed
sim.hospital=2; % fitting to both incidence and hosp admissions
sim.hospitalReportingLag=0;

%% Load population data
sim.populationHandle = fullfile(sim.dataFolder, 'population2020.txt');
populationTable = readtable(sim.populationHandle);
sim.population = populationTable.Population';

%% Load k-value
sim.kValueHandle = fullfile(sim.dataFolder, 'kvalue.csv');
kvalTable = readtable(sim.kValueHandle);
sim.kval = kvalTable.Var1;
if length(sim.kval) ~= length(sim.dateList)
    %disp("Error: kValue List not equal length of date list")
    if length(sim.kval)>length(sim.dateList)
        sim.kval=sim.kval(1:length(sim.dateList));
    else
        for j=length(sim.kval)+1:length(sim.dateList)
            sim.kval=[sim.kval;1];
        end
    end
end

%% Load demographics
% fullFileName = fullfile(pwd,sim.folder,'data/demographic.txt');
% ds = loadMultiVarFormat(fullFileName);
% sim.q = ds.q;
% sim.mu = ds.mu;

sim.q = 0; % Set quarantine parameter, modify as needed.
sim.mu = zeros(1,16); % Set mortality rates, modify as needed.

%% Load pax data
% sim.paxHandle = fullfile(pwd, sim.folder, 'data/pax.txt');
% paxTable = readtable(sim.paxHandle);
% paxIndex = ismember(datenum(datetime(paxTable.date)), sim.dateNumList);
% sim.pax = paxTable.pax(paxIndex);

%% Load omega data
sim.omegaHandle = fullfile(sim.dataFolder, 'omega.txt');
omegaTable = readtable(sim.omegaHandle);
omegaDateNum = datenum(datetime(omegaTable.date));
omegaIndex = ismember(omegaDateNum,sim.dateNumList);
fieldNames = {'omega', 'omega2', 'omega3', 'omega4'};
for i = 1:length(fieldNames)
    sim.(fieldNames{i}) = omegaTable.omega(omegaIndex);
end

%% Load contact modifiers
contactModifersTable = readtable(fullfile(sim.dataFolder, 'contactModifiers_weekly.txt'));
sim.contactModifers = contactModifersTable{:,2:4};

%% Load contact matrices
fileNames = {'italy_contactHome.txt', 'italy_contactWork.txt', 'italy_contactOtherPlaces.txt', 'italy_contactSchool.txt'};
fieldNames = {'H', 'W', 'O', 'S'};
for i = 1:length(fileNames)
    fullFileName = fullfile(sim.dataFolder, fileNames{i});
    sim.(fieldNames{i}) = load(fullFileName);
end

%% Load perturbation matrices
prefixes = {'SF', 'WF', 'OF'};
prefixNumbers = [3, 2, 1];
for i = 1:length(prefixes)
    for j = 1:length(prefixNumbers)
        fileName = ['contact', prefixes{i}, num2str(prefixNumbers(j)), '.txt'];
        fullFileName = fullfile(sim.dataFolder, fileName);
        sim.(sprintf('%s%d', prefixes{i}, prefixNumbers(j))) = load(fullFileName);
    end
end

%% Fit proportions
%sim.fitProportions_inc = [0 0.05 0.6 1];
sim.fitProportionsLabel_inc = 'Variable';
sim.fitProportionsHandle_inc = fullfile(sim.dataFolder, 'fitPropInc_Var_July2020.csv');
fitPropTable= readtable(sim.fitProportionsHandle_inc);
sim.fitProportions_inc = table2array(fitPropTable(:,2:end));

% sim.fitProportions_hosp = [0 0 0 0.4];
% sim.fitProportionsLabel_hosp = '0-0-0-.4';
sim.fitProportionsLabel_hosp = 'Variable';
sim.fitProportionsHandle_hosp = fullfile(sim.dataFolder, 'fitPropHosp_Var_July2020.csv');
fitPropTable= readtable(sim.fitProportionsHandle_hosp);
sim.fitProportions_hosp = table2array(fitPropTable(:,2:end));

%% Comorbidities
comorb_handle=fullfile(sim.dataFolder, 'p1_Italy.txt');
comorbTable = readtable(comorb_handle);
sim.p1=table2array(comorbTable);


%% Load vaccination data
sim.vaccinateHandle = fullfile(sim.dataFolder, 'vaccinate.txt');
vaccinateTable = readtable(sim.vaccinateHandle);
sim.vaccinate = vaccinateTable.vaccinate;

%% Load vaccination coverage
sim.vaccineCoverageHandle = fullfile(sim.dataFolder, 'vaccineCoverage.txt');
sim.vaccineCoverageLabel=0; % Coverage data given as absolute numbers, not percentages
sim.vacTargetPopulationHandle = fullfile(sim.dataFolder, 'vaccineTargetPopulation.txt');
vacTargetTable = readtable(sim.vacTargetPopulationHandle);
vacTargetData = table2array(vacTargetTable(:,2:end));
vacTargetDates = datenum(vacTargetTable.Var1);
sim.vacTargetPopulation = vacTargetData(ismember(vacTargetDates,sim.dateNumList),:);
sim.newVacTarget = vacTargetData(ismember(vacTargetDates,sim.dateNumList),:);

%% Set handles to vaccination rates (here we don't calculate sigmas from vaccination rate data, but import them from previously-calculated files)
sim.sigmaLabel="Import";
sim.sigmasHandles = [ ...
    string(fullfile(sim.dataFolder, 'sigma1.txt')), ...
    string(fullfile(sim.dataFolder, 'sigma2.txt')), ...
    string(fullfile(sim.dataFolder, 'sigmaB.txt')) ...
];

sim.sigma1=sim.get_sigma(1);
sim.sigma2=sim.get_sigma(2);
sim.sigmaB=sim.get_sigma(3);


%% setting region, deltaWave, omicron, etc: 
sim.region = 'Italy';
sim.alphaWaveStart= 738145; %First reported case 20-Dec-2020 in Italy https://tg24.sky.it/mondo/2020/12/20/covid-variante-inglese
sim.deltaWaveStart=738284; % First reported case 08-May-2021  in Italy https://ilsicilia.it/covid-oltre-30-casi-variante-delta-in-sicilia-2-le-zone-rosse/
sim.omicronWaveStart = 738487; % First reported case 27-Nov-2021 in Italy https://www.ilgiornale.it/news/cronache/chi-paziente-zero-variante-omicron-1992281.html
sim.acuteCareCap = 156216; % Acute care bed capacity prior to COVID-19 in Italy https://pmc.ncbi.nlm.nih.gov/articles/PMC8632742/


%ICU_capacity = 8550; %ICU bed surge capacity in Italy during 1st wave 




% sim.deltaWave = 738322; %15-Jun-2021
% sim.omicronWave = 738502; %12-Dec-2021

%% Fit properties
%% Fit properties
sim.curvaturePenaltyPercent = 1e-3; % Sets the penalty percentage for the curvature in the fit.
sim.curvaturePenalty = 1e1;         % Sets the absolute value of the penalty for the curvature.
sim.populationFactor = 6;            % Sets the scaling factor for the population used in the simulation.
sim.timeLimit = 40;                  % Sets the time limit per fit window (minutes).
sim.windowLoopStart = 1;           % Sets the starting index for the sliding fit window.
sim.windowLoopEnd = 190;             % Sets the ending index for the sliding fit window.
sim.windowSize = 16;                 % Sets the size of the sliding fit window.
sim.windowStep = 4;                  % Sets the step size for moving the sliding fit window.
sim.upperBound = 12;                  % Sets the upper bound for the optimization algorithm.
sim.lowerBound = 0.1;                % Sets the lower bound for the optimization algorithm.
sim.maxStallTime = 60*5;             % Sets the maximum time to allow the optimizer to stall (in seconds).
sim.maxGenerations = 240;            % Sets the maximum number of generations for the optimization algorithm.


%% file handles
sim.forceFitHandle = fullfile(sim.dataFolder, 'forceFit.txt');
sim.exitFlagHandle = fullfile(sim.dataFolder, 'exitFlags.txt');

%% Set booster dose decay
sim.boosterDeacyFlag = false;

%%
sim.beta3=0.0694;

%%
sim.LOS_Q1=8;
sim.LOS_median=12;
sim.LOS_Q3=21;

% sim.LOS_Q1=9;
% sim.LOS_median=14;
% sim.LOS_Q3=24;

%% Run simulation
sim.simulate();

end
