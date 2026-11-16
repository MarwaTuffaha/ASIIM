% Add the shared ASIIM functions to the MATLAB path
italyFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(italyFolder, '..', 'functions'));

clc

% Define simulation and load Italy data
sim = ItalySetup();
plotBothData(sim)

%% Plot fit with the initial k-values
sim.fitPlot

%% Choose fitting window (same settings as the original example)
sim.windowLoopStart = 50;
sim.windowLoopEnd = 102;
% sim.fitPlot(true,[sim.windowLoopStart,sim.windowLoopEnd])

%% Fit properties (unchanged from the original example)
sim.timeLimit = 3;             % minutes per window
sim.windowSize = 16;
sim.windowStep = 2;
sim.upperBound = 1.5;
sim.lowerBound = 0.5;
sim.monitorPlot = false;
sim.maxStallTime = 60*10;
sim.useParallel = true;
sim.maxGenerations = 300;
sim.curvaturePenaltyPercent = 1e-1;

sim.fit

%% Simulate and plot the fitted result
sim.simulate
sim.fitPlot

%% Optional results (require the corresponding shared plotting functions)
sim.immunity
sim.compartments
