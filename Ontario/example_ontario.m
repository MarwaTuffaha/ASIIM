% Add the shared ASIIM functions to the MATLAB path
ontarioFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(ontarioFolder, '..', 'functions'));


% Define simulation
% Loading the data into the object
sim = ontarioSetup();
plotBothData(sim)

%% Plot Fit
sim.fitPlot

%% Showing how the plot looks like when specifying a fitWindow

sim.windowLoopStart = 2;  %choose first window that has not been fit to data
sim.windowLoopEnd = 107;  %choose last window with data

%sim.fitPlot(true,[sim.windowLoopStart,sim.windowLoopEnd])

% Fit properties
        sim.timeLimit = 0.8;       %per loop window fit time limit in minutes
        sim.windowSize = 16;    %looped fit window size
        sim.windowStep = 2;     %looped fit window step size
        sim.upperBound = 4;       %Upper Bound on the k-value (bigger = slower)
        sim.lowerBound = 1e-2;    %soultion time may be large when lb ~ 0
        sim.monitorPlot = false;  %Update fitPlot every time a new best k-value is found
        sim.maxStallTime = 60*10;    %The algorithm stops if there is no improvement in the best fitness value (seconds)
        sim.useParallel = true;   %compute GA algorithm in parallel , test!
        sim.maxGenerations = 300; %MaxGenerations option determines the maximum number of generations the genetic algorithm takes
        sim.curvaturePenaltyPercent=1e-1;

% Run fit
sim.kValueHandle = fullfile(sim.folder, 'results', 'kvalue_calibrated.csv'); % Preserve input kvalue.csv
sim.fit


%%
% Simulate
sim.simulate

% Plot Fit
sim.fitPlot

%% Immunity Results
sim.immunity  
