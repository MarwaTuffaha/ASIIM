classdef asiim < matlab.mixin.Copyable
    % ASSIM The experimental class for Covid-Assim
    %
    % This class implements the Covid-Assim simulation model, which is used
    % to analyze and predict the spread of COVID-19 within a population.
    % The class includes properties for configuring the simulation, storing
    % data, and managing model parameters, as well as methods for running
    % simulations, fitting models to data, and plotting results.
    %
    % Properties:
    %   simNum              - Simulation number (default: 1)
    %   experimentDir       - The full path of the experinment root folder
    %   experimentName      - Name of the experiment (default: 'default')
    %   folder              - Folder where the experiment data is stored
    %   description         - Description of the simulation
    %   shortTitle          - Short title for the experiment
    %   threshold           - Threshold to set compartment to zero if value < threshold
    %   saveDiagnosticFlag  - Flag to save diagnostic data (default: true)
    %   diagnosticFolder    - Folder for diagnostic data (default: '/fitDiagnostic/')
    %   calenderStopDate    - Calendar stop date for the simulation
    %   dateListHandle      - File location of date list in yyyy-mm-dd format
    %   dateList            - List of dates
    %   dataFolder          - Folder for data (default: '/data/')
    %   hospitalFile        - File name of hospitalization data (default: 'hospitalizationData.csv')
    %   dataTypeLabel_inc   - Flag indicating if incidence data type is prevalence (2), cummulative (1) or daily (0) (default: 0)
    %   dataTypeLabel_hosp  - Flag indicating if hospitalization data type is prevalence (2), cummulative (1) or daily (0) (default: 0)
    %   hospDates           - Dates of hospital data
    %   hospCases           - Cases of hospital data
    %   hospCasesSmoothLabel- Label indicating whether the hosp data is already smoothed (1) or not (0)
    %   hospitalReportingLag - Number of days to shift hospital data for fit (default: 3)
    %   incidenceFile       - File name of incidence data (default: '')
    %   incidenceDates      - Dates of incidence data
    %   incidenceCases      - Cases of incidence data
    %   country             - Country of the simulation
    %   region              - Region of the simulation (default: 'Ontario')
    %   populationHandle    - Location of population file
    %   population          - Population data
    %   beta3               - Beta parameter for the model (default: 0.08)
    %   hospital            - Flag indicating if hospital data is used (default: true)    0=incidence only | 1=hospital only  |  2=Both data
    %   experimentDescription - Description of the experiment
    %   H                   - Home contact matrix
    %   S                   - School contact matrix
    %   W                   - Work contact matrix
    %   O                   - Other contact matrix
    %   SF1, SF2, SF3, SF4  - Contact matrices for different settings
    %   WF1, WF2, WF3, WF4  - Contact matrices for different settings
    %   OF1, OF2, OF3, OF4, OF5, OF6, OF7 - Contact matrices for different settings
    %   deltaBool           - Flag for delta parameter (default: 0)
    %   agingBool           - Flag for aging parameter (default: 0)
    %   ageGroupInfected    - Age groups to infect in initial conditions
    %   IO                  - Initial infected (default: 1/10)
    %   SIC                 - Initial conditions for susceptible-infectious compartments
    %   mu                  - Mortality rates
    %   q                   - Quarantine parameter
    %   omegaHandle         - Handle for omega parameter
    %   omegaDate           - Dates for omega parameter
    %   omega               - Omega parameter values
    %   omegaElement        - Omega parameter elements
    %   omega4, omega3, omega2 - Omega parameter values
    %   omegaLabel          - Label for omega parameter
    %   paxHandle           - Handle for Paxlovid data
    %   pax                 - Time-varying protection rates from Paxlovid (daily)
    %   deltaCoverageLabel  - Label for delta coverage (default: 'Base')
    %   vaccineCoverageData - Vaccine coverage data (default: 0)
    %   vaccineCoverageHandle - Handle for vaccine coverage data
    %   vaccineCoverageLabel - Label for whether coverage data is given as absolute numbers (0) or percentage (1)
    %   maxcoverage         - Maximum vaccine coverage
    %   maxcoverageTwo      - Second maximum vaccine coverage
    %   maxcoverageBoost    - Maximum booster vaccine coverage
    %   vaccinateData       - Vaccination data (default: 0)
    %   vaccinateHandle     - Handle for vaccination data
    %   vaccinate           - Vaccination data
    %   vaccinateLabel      - Label for vaccination data (default: 'vaccination')
    %   boosterDeacyFlag    - Flag for booster decay (default: false)
    %   doseLabel           - Label for dose
    %   vacTargetPopulationHandle - Handle for vaccine target population data
    %   vacTargetPopulation - Vaccine target population data
    %   newVacTarget        - New vaccine target population data
    %   vacTargetPopulationLabel - Label for vaccine target population (default: 1)
    %   vacTargetPopulationLabelVerbose - Verbose label for vaccine target population (default: 'HCW')
    %   contactModifers     - Contact modifiers data
    %   contactLabel        - Label for contact data
    %   contactDateLabel    - Label for contact date data
    %   rho                 - Vaccine efficacy (default: 0.9)
    %   rhoLabel            - Label for rho (default: '90%')
    %   epsilon             - Vaccine protection from single dose (default: 0.7)
    %   epsilonLabel        - Label for epsilon (default: '70%')
    %   qq                  - Vaccine protection from severe disease (default: 0.92)
    %   qLabel              - Label for qq (default: '3')
    %   qq_SLabel           - Label for whether we apply the protection from severe disease to susceptible classes (1) or only to V1 compartments ( default: 0)
    %   vaccineTypeLabel    - Label for vaccine type (default: '3')
    %   gamma2              - Recovery rate for compartment 2 (default: 1/5)
    %   gamma3              - Recovery rate for compartment 3 (default: 1/10)
    %   gamma4              - Recovery rate for compartment 4 (default: 1/15)
    %   delta4              - Disease-induced mortality for compartment 4 (default: 0)
    %   delta3              - Disease-induced mortality for compartment 3 (default: 0)
    %   alpha               - Alpha parameters
    %   alphaLabel          - Label for alpha (default: '[1 .67 .33]')
    %   forceFitBool        - Flag for force fit (default: true)
    %   forceFitData        - Data for force fit (default: 0)
    %   kValueHandle        - Handle for k-value data
    %   kval                - K-value data
    %   kvalLabel           - Label for k-value (default: {'Force Fit'})
    %   fitProportionsHandle_inc    - Handle for fit proportions data for incidence data
    %   fitProportions_inc          - Fit proportions data for incidence data
    %   fitProportionsLabel_inc     - Label for fit proportions for incidence data
    %   fitProportionsHandle_hosp   - Handle for fit proportions data for hospitalizations data
    %   fitProportions_hosp         - Fit proportions data for hospitalizations data
    %   fitProportionsLabel_hosp    - Label for fit proportions for hospitalizations data   
    %   forceFitHandle      - Handle for force fit data
    %   voc                 - VOC parameter (default: 1.2)
    %   c                   - Additional parameter
    %   A                   - Additional parameter
    %   mitagationWindowLengthElement - Length of mitigation window element
    %   time                - Simulation time
    %   yout                - Simulation output
    %   timeDetail          - Detailed simulation time
    %   youtDetail          - Detailed simulation output
    %   youtUB              - Upper bound estimate for yout
    %   youtLB              - Lower bound estimate for yout
    %   yint                - Custom starting initial conditions
    %   useInitialConditions - Flag to use initial conditions (default: false)
    %   yintDetail          - Initial conditions for seroprevalence equations
    %   useDetailInitialConditions - Flag to use detailed initial conditions (default: false)
    %   fitDiagnositcPlot   - Flag to update diagnostic plot during fitting (default: false)
    %   populationFactor    - Population factor (default: 3)
    %   timeLimit           - Time limit per loop window fit (default: 60 minutes)
    %   windowLoopStart     - First window that has not been fit to data (default: 1)
    %   windowLoopEnd       - Last window with data (default: 105)
    %   windowSize          - Looped fit window size (default: 17)
    %   windowStep          - Looped fit window step size (default: 2)
    %   upperBound          - Upper bound on the k-value (default: 6)
    %   lowerBound          - Lower bound on the k-value (default: 1e-2)
    %   monitorPlot         - Flag to update fit plot when a new best k-value is found (default: false)
    %   maxStallTime        - Maximum stall time in seconds (default: 600)
    %   useParallel         - Flag to use parallel computation (default: true)
    %   maxGenerations      - Maximum number of generations for GA algorithm (default: 300)
    %   curvaturePenalty    - Initial curvature penalty (default: 1e-1)
    %   curvaturePenaltyPercent - Curvature penalty as percent of objective function (default: 1e-2)
    %   fitVariance         - Flag to fit variance (default: false)
    %   variancePenaltyPercent - Variance penalty as
    %   cut_Omicron         - Flag to indicate whether to cut data after Omicron wave (1) or not (default:0)
    %   morePoints          - Number of data points to keep after the Omicron wave data to improve the fit at the last few points
    %   S_fLabel            - Label indicating whether we scale the force of infection by the susceptible-saturation factor, 0=off (default), 1=on 
    %   alphaWaveStart      - First reported alpha case
    %   deltaWaveStart      - First reported delta case
    %   omicronWaveStart    - First reported omicron case
    %   acuteCareCap        - Acute care capacity
    %   LOS_Q1              - 25th percentile length of stay in the hospital
    %   LOS_median          - Median length of stay in the hospital
    %   LOS_Q3              - 75th percentile length of stay in the hospital
    %   LOS_dist            - distribution to use for the hospital length of stay: PERT or lognormal
    %   alphaInput          - susceptibility parameters



    properties
        simNum = 1;
        experimentDir = '';
        experimentName string = 'default'
        folder string = '/experiments/'; %folder where the experiment data is stored
        description string = 'simulation description'
        shortTitle string = 'shortTitle';
        cut_Omicron=0;
        morePoints=0;
        S_fLabel=0;
        
        threshold = 10e-8; % Set compartment to zero if value < threshold

        %% fit diagnostic properties
        saveDiagnosticFlag = true;
        diagnosticFolder = '/fitDiagnostic/'

        %% timming properties
        calenderStopDate = '2020/1/1';

        dateListHandle %dateList File location yyyy-mm-dd
        dateList(:,3)

        %% Data
        dataFolder string = '/data/'
        hospitalFile = 'hospitalizationData.csv';
        dataTypeLabel_inc = 0;
        dataTypeLabel_hosp = 0;

        hospDates(:,1) double
        hospCases(:,1) double
        hospCasesSmoothLabel=0;
        hospitalReportingLag = 3;%Number of days to shift hospital data for fit

        incidenceFile = ''
        incidenceDates(:,1) 
        incidenceCases(:,1) 

        country string
        region string = 'Ontario';
        %totalPopulation(1,1) int16

        %% Simulation properties

        populationHandle %{mustBeFile} Location of population file
        population(1,16)

        beta3 = 0.0545;

        hospital; %use hospital data (this is requierd atm)
        experimentDescription string

        H(16,16) double %Home Contact Matrix
        S(16,16) %School Contact Matrix
        W(16,16) %Work Contact Matrix
        O(16,16) %Other Contact Matrix
        SF1(16,16)
        SF2(16,16)
        SF3(16,16)
        SF4(16,16)
        WF1(16,16)
        WF2(16,16)
        WF3(16,16)
        WF4(16,16)
        OF1(16,16)
        OF2(16,16)
        OF3(16,16)
        OF4(16,16)
        OF5(16,16)
        OF6(16,16)
        OF7(16,16)
        deltaBool =0;
        agingBool =0;
        ageGroupInfected(1,:) = [5,6,7,8,9,10];%age groups to infect in initial conditions
        IO(1,1) = 1/6; %Initial infected
        SIC = ones(1,16);

        mu(1,16) = 0;
        q(1,1) = 0;

        omegaHandle
        omegaDate
        omega double
        omegaElement
        omega4(:,1)
        omega3(:,1)
        omega2(:,1)
        omegaLabel = 'time variable';

        paxHandle
        pax(:, 1) %contains the time-varying protection rates from Paxlovid (daily)

        deltaCoverageLabel = 'Base';

        vaccineCoverageData = 0;
        vaccineCoverageLabel
        vaccineCoverageHandle
        maxcoverage(1,:)
        maxcoverageTwo(1,:)
        maxcoverageBoost(1,:)

        vaccinateData = 0;
        vaccinateHandle
        vaccinate(:,1)
        vaccinateLabel = 'vaccination';

        boosterDeacyFlag = false;

        doseLabel = 0;

        vacTargetPopulationHandle
        vacTargetPopulation(:,16)
        newVacTarget(:,16)
        vacTargetPopulationLabel = 1;
        vacTargetPopulationLabelVerbose = 'HCW';

        contactModifers(:,3)
        contactLabel
        contactDateLabel

        rho(1,1) = 0.9; %vaccine efficacy
        rhoLabel = '90%' %was 90 for CA, use 80 for Hun

        epsilon = .7 %vaccine protection from singe dose
        epsilonLabel = '70%';

        qq = 0.92 %vaccine protection from severe disease
        qLabel = '3';

        vaccineTypeLabel ='3';

        % recovery (per day)
        gamma2 = 1/5;
        gamma3 = 1/10;
        gamma4 = 1/15;

        % disease-induced mortality (=0) %turned Off
        delta4 = 0;
        delta3 = 0;

        alphaLabel ='[1 .67 .33]';

        forceFitBool logical = 1;
        forceFitData = 0;

        kValueHandle
        %kval
        kvalLabel = {'Force Fit'};

        fitProportionsHandle_inc;
        fitProportions_inc;%(1,4) = [0 0 0 .2];
        fitProportionsLabel_inc

        fitProportionsHandle_hosp;
        fitProportions_hosp;%(1,4) = [0 0 0 .2];
        fitProportionsLabel_hosp        

        forceFitHandle

        voc = 1.2;

        kval(:,1)

        boosterDecayFlag = false;

        c
        A

        mitagationWindowLengthElement

        %ODE results
        time
        yout
        timeDetail
        youtDetail
        youtUB %Upper bound estimante for yout
        youtLB %Lowwer bound estimate for yout

        yint %custom starting inital condidtions
        useInitialConditions = false;

        yintDetail %initial conditions for the seroprevelence equations
        useDetailInitialConditions = false;

        deltaWave = 738335+7; % Maybe end of wave 05-Jul-2021
        omicronWave = 738505; % midpoint of wave

        alphaWaveStart= 738126; %First reported case 01-Dec-2020 in Ontatio
        deltaWaveStart=738224; % First reported case 09-Mar-2021  in Ontario
        omicronWaveStart = 738490; % First reported case 30-Nov-2021 in Ontario
        acuteCareCap= 29743+2300; % Acute Care capacity in Ontario https://link.springer.com/article/10.1007/s11538-023-01130-x

        qq_SLabel=0;

        % These are the values for Ontario hospital LOS https://www.cmaj.ca/content/193/12/E410
        LOS_Q1=3.6;
        LOS_Q3=18.9;
        LOS_median=8.7;
        LOS_dist="lognormal";

    end

    %% Fitting properties
    properties
        fitDiagnositcPlot = false; %Update Diagnositc plot during fitting

        populationFactor = 3;
        timeLimit = 60;       %per loop window fit time limit in miniutes
        windowLoopStart = 1;  %choose first window that has not been fit to data
        windowLoopEnd = 50;  %choose last window with data
        windowSize = 17;    %looped fit window size
        windowStep = 2;     %looped fit window step size
        upperBound = 6;       %Upper Bound on the k-value (bigger = slower)
        lowerBound = 1e-2;    %soultion time may be large when lb ~ 0
        monitorPlot = false;  %Update fitPlot every time a new best k-value is found
        maxStallTime = 60*10;    %The algorithm stops if there is no improvement in the best fitness value (seconds)
        useParallel = true;   %compute GA algorithm in parallel , test!
        maxGenerations = 300; %MaxGenerations option determines the maximum number of generations the genetic algorithm takes

        curvaturePenalty = 1e-1; %Initial curvature penalty
        %curvature penalty (i.e. Regularization Hyperparameter)
        %With the assumption that the public changes in the level of mitagation
        %is 'sticky' in the economic sense. Rapid changes to the fit k-vlaue are
        %penalized using an approximation of curvature of the k-value vector.
        %Choose a value that is small relative to the error.
        %NOTE: this is also reset after obj is normalized.
        curvaturePenaltyPercent = 1e-2; %curvature penalty as percent of objective function
        fitVariance = false;
        variancePenaltyPercent = 1e-1;

        fitTime = 0;
        objectiveFunctionValue = Inf;
        exitFlagHandle

        p1 = [0.9792    0.0195    0.0013
            0.9703    0.0280    0.0017
            0.9714    0.0270    0.0017
            0.9638    0.0338    0.0024
            0.9474    0.0486    0.0040
            0.9237    0.0693    0.0070
            0.8972    0.0915    0.0113
            0.8656    0.1168    0.0176
            0.8070    0.1625    0.0305
            0.7566    0.1979    0.0456
            0.6902    0.2415    0.0683
            0.6002    0.2963    0.1035
            0.5032    0.3468    0.1500
            0.4095    0.3840    0.2065
            0.3247    0.4040    0.2713
            0.2152    0.3738    0.4111]'; %Canada coMorbidity (infection fraction from S1 to I2,I3,I4)
        
        sigmaLabel="Calculate";  % Label that specifies whether the vaccination rates are calculated from vaccine coverage data ("Calculate") or imported ("Import")
        sigmasHandles       % Handles to vaccination rates files (if to be imported)
        sigma1,sigma2,sigmaB              % Vaccination rates (if to be imported)

        alphaInput(1,3) double = [1 2/3 1/3];
        

    end
    properties (Constant)
        initPopSize = 15; %Changes might not work
       
        age = 2.5:5:80; % Mid point of age interval  (in years)
        age_interval = [5*ones(1,15) 25]; % age intervals (in years)



        % US population data per age class used for perFemales in fecundity
        % caculation
        USmales = [10319427, 10389638, 10579862, 11303666,...
            11014176, 10635591, 9996500, 10042022,...
            10393977, 11209085, 10933274, 9523648,...
            8077500, 5852547, 4243972, 3182388+2294374+1273867+424387+82263+9162];

        USfemales = [9881935,9959019,10097332, 10736677,...
            10571823,10466258,9965599,10137620,...
            10496987,11499506,11364851,10141157,...
            8740424,6582716,5034194,4135407+3448953+2346592+1023979+288981+44202];

        % fecundity (originally per year scaled to per day)
        fRaw = [0, 0, 0, 1.2150e-04,...
            1.0328e-02,...
            3.8922e-02,...
            5.3160e-02,...
            5.2850e-02,...
            2.4996e-02,...
            5.1066e-03,...
            3.3678e-04,...
            3.0843e-05, zeros(1,4)]/365;

        % transition from E to I
        kappaZero = 1.5; % average duration in an exposed class (days)
        kappa1 = 0; % zeroing out the rate at index 1

        I0 = 0;
        t3 = 200;

        f = 0; %demographics are off
        

    end
    properties (Hidden)
        daysbetweendoses = 112;%Depricated
        outLoop
        sim
        fitStop = 5*365;
        deltaCoverageBoost
        sigmaError
        short
    end
    properties (Dependent)
        totalPopulation(1,1)
        deltaCoverage(1,:)
        deltaCoverageTwo
        %
        forceFit(:,1)
        dateNumList(:,1)
        mitagationWindowLength(:,1)
        mitagationWindowLengthShort(:,1)
        fitStopVec
        vaccinationStart
        vaccinationBoostStart
        vaccineCoverage
        fitCurvature
        error
        mse
        kvalVar
        tmax
        N
        UStotal
        PerFemales
        alpha(1,3)
        %         f
        gamma
        beta
        p2 % infection fraction from S2 to I2,I3,I4
        p3 % infection fraction from S3 to I2,I3,I4
        delta
        kappa
        tau
        normalTest
        pastPredictionDate
        stopDate

        simulationLogFilePath
        errorLogFilePath
        

    end

    methods

        function obj = asiim()
            % Constructor to initialize the object

            % Check if results folder exists
            if ~exist(obj.folder, 'dir')
                % Create results folder if it doesn't exist
                mkdir(obj.folder);
                disp(['Created folder: ', obj.folder]);
            end

%             obj.N = N;
%             obj.equations = equations;
%             obj.Index = struct(); % Initialize Index as an empty struct
%             
%             % Initialize compartmentSize property
%             obj.compartmentSize = struct();
%             obj.compartmentSize.original = [8, 8, 8, 8, 8, 5, 5, 5];
%             obj.compartmentSize.detailed = [12, 8, 8, 8, 8, 5, 5, 7];
        end

        %% Log files
        function path = get.simulationLogFilePath(obj)
             if ~exist(fullfile(obj.folder,'/logFiles'), 'dir')
                % Create results folder if it doesn't exist
                mkdir(fullfile(obj.folder,'/logFiles'));
                disp(['Created folder: ', obj.folder]);
            end
            path = fullfile(obj.folder, 'logFiles/simulationLog.txt');
        end

        function path = get.errorLogFilePath(obj)
             if ~exist(fullfile(obj.folder,'/logFiles'), 'dir')
                % Create results folder if it doesn't exist
                mkdir(fullfile(obj.folder,'/logFiles'));
                disp(['Created folder: ', fullfile(obj.folder,'/logFiles')]);
            end
            path = fullfile(obj.folder,'/logFiles/errorLog.txt');
        end

        %% error checks - for text files
        function checkData(obj)
            dateError = dateListErrorCheckFunc(obj);

            caseError = 0;
            if obj.hospital
                caseError = checkCaseData(obj.hospCases,obj);
            end

            if dateError||caseError
                disp('Errors found: check log files')
            end
        end

        %% properties and validation
        function totalPopulation = get.totalPopulation(obj)
            totalPopulation = sum(obj.population);
        end
        function deltaCoverage = get.deltaCoverage(obj)
            deltaCoverage=diff(obj.vaccineCoverage(:,1));
            if ~isempty(deltaCoverage)
                deltaCoverage(end+1)=deltaCoverage(end);
            end
        end
        function deltaCoverageTwo = get.deltaCoverageTwo(obj)
            deltaCoverageTwo=diff(obj.vaccineCoverage(:,2));
            if ~isempty(deltaCoverageTwo)
                deltaCoverageTwo(end+1)=deltaCoverageTwo(end);
            end
        end

        %this should use number of mitagation windows?
        function set.windowLoopEnd(obj, value)
            obj.windowLoopEnd = obj.validateWindowLoopEnd(value);
        end
        function ret = validateWindowLoopEnd(obj,value)
            len = length(obj.dateNumList);
            if obj.windowLoopEnd >= len
                fprintf('Invalid value for windowLoopEnd: %d. Set to one less then length of dateNumList: %d.\n', ...
                    obj.windowLoopEnd, len-1);
                ret = len-1;
            else
                ret = value;
            end
        end

        % %         function deltaCoverageBoost = get.deltaCoverageBoost(obj)
        % %             deltaCoverageBoost=diff(obj.vaccineCoverage(:,3));
        % %             deltaCoverageBoost(1:end+1)=[0;deltaCoverageBoost];
        % %         end
        function value = get.stopDate(obj)
            value = datenum(obj.calenderStopDate);
        end

        function set.fitTime(obj,value)
            writematrix(value,fullfile(obj.experimentDir(),obj.folder,'/results/totalfittime.txt'))
            obj.fitTime = value;
        end
        % %         function dateList = set.dateListHandle(obj,value)
        % %             dateTable = readtable(obj.dateListHandle);
        % %             dateListDT = datetime(dateTable.Var1,'InputFormat','yyyy-MM-dd');
        % %             dateList = [dateListDT.Day dateListDT.Month dateListDT.Year];
        % %         end
        function forceFit = get.forceFit(obj)
            if length(obj.forceFitData) == 1
                fitTable = readtable(obj.forceFitHandle);
                forceFit = fitTable.Var2;
            else
                forceFit = obj.forceFitData;
            end
            if obj.cut_Omicron==1
                forceFit=forceFit(1:length(obj.dateList));
            end
            if length(forceFit) ~= length(obj.dateList)
                disp("Error: forceFit list not equal length of date list")
            end
        end
        function dateNumList = get.dateNumList(obj)
            dateNumList = datenum(obj.dateList(:,[3,2,1]));
        end
        function mitagationWindowLength = get.mitagationWindowLength(obj)
            mitagationWindowLength = diff(obj.dateNumList);
        end
        function mitagationWindowLengthShort = get.mitagationWindowLengthShort(obj)
            mitagationWindowLengthShort = diff(obj.dateNumList(obj.vaccinate==1));
        end
        function fitStopVec = get.fitStop(Ont)
            fitStopVec = ...
                1:min(Ont.windowSize,Ont.windowStep):...
                (Ont.dateNumList(Ont.windowLoopEnd)-Ont.dateNumList(Ont.windowLoopStart));
        end
        function vaccinationStart = get.vaccinationStart(obj)
            vaccinationStart = obj.dateNumList(find(obj.vaccinate, 1, 'first'));
        end
        function vaccinationBoostStart = get.vaccinationBoostStart(obj)
            shortList = obj.dateNumList(logical(obj.vaccinate));
            vaccinationBoostStart = shortList(find(obj.vaccineCoverage(:,3)~=0, 1, 'first'))-7;%one window early?
        end

        function vaccinate = get.vaccinate(obj)

            % Check the length of the vaccinateData property
            if length(obj.vaccinateData) == 1
                % If vaccinateData has only one element, we use the file handle.
                % Read the vaccination data from the file specified by the vaccinateHandle property.
                % vaccinateTable is a table that stores the loaded vaccination data.
                vaccinateTable = readtable(obj.vaccinateHandle);
                
                % Convert 'vaccinate' column to logical array
                vaccinate = logical(vaccinateTable.vaccinate);
            else
                % If vaccinateData contains more than one element, it is assumed to be the actual data.
                % In this case, directly use the vaccinateData as the vaccination table.
                vaccinate = logical(obj.vaccinateData);
            end



            % Find indices where date numbers exceed the last incidence or hospitalization date
            if obj.hospital==1
                index = find(obj.dateNumList > obj.hospDates(end));
            else
                index = find(obj.dateNumList > obj.incidenceDates(end));
            end

            % Set 'vaccinate' to false for these indices
            vaccinate(index) = false(length(vaccinate(index)), 1);
            if obj.cut_Omicron==1
                vaccinate=vaccinate(1:length(obj.dateNumList));
            end
            % Check if the length of 'vaccinate' array matches with 'dateNumList'
                        
            if (length(vaccinate) - length(obj.dateNumList)) ~= 0
                % Display length discrepancy if any
                disp(['vaccinate list length error by: ' ...
                    num2str(length(vaccinate) - length(obj.dateNumList))])

            end
%         function vaccinate = get.vaccinate(obj)
% 
%             % Check the length of the vaccinateData property
%             if length(obj.vaccinateData) == 1
%                 % If vaccinateData has only one element, we use the file handle.
%                 % Read the vaccination data from the file specified by the vaccinateHandle property.
%                 % vaccinateTable is a table that stores the loaded vaccination data.
%                 vaccinateTable = readtable(obj.vaccinateHandle);
%                 
%                 % Convert 'vaccinate' column to logical array
%                 vaccinate = logical(vaccinateTable.vaccinate);
%             else
%                 % If vaccinateData contains more than one element, it is assumed to be the actual data.
%                 % In this case, directly use the vaccinateData as the vaccination table.
%                 vaccinate = logical(obj.vaccinateData);
%             end
% 
% 
% 
%             % Find indices where date numbers exceed the last hospitalization date
%             index = find(obj.dateNumList > obj.hospDates(end));
% 
%             % Set 'vaccinate' to false for these indices
%             vaccinate(index) = false(length(vaccinate(index)), 1);
% 
%             % Check if the length of 'vaccinate' array matches with 'dateNumList'
%             if (length(vaccinate) - length(obj.dateNumList)) ~= 0
%                 % Display length discrepancy if any
%                 disp(['vaccinate list length error by: ' ...
%                     num2str(length(vaccinate) - length(obj.dateNumList))])
%             end
        end

        function sigma = get_sigma(obj,dose)
            % ---- Inputs ----
            sigmaDailyTable = readtable(obj.sigmasHandles(dose));  % has sigmaDailyTable.date + age cols
            [sigmaDailyTable, isError] = uniqueVaccineDatesFunc(sigmaDailyTable, obj.dataFolder);
            if isError
                logAndThrowError(obj, 'Daily sigma dates are not unique.');
            end
            
            weekDates = obj.dateNumList(logical(obj.vaccinate));   % datenum weekly list (your calibration list)
            
            % ---- Convert daily table to arrays ----
            dDaily = datenum(sigmaDailyTable.Var1);
            Xdaily = table2array(sigmaDailyTable(:,2:end));  % [Ndaily x Nage]
            
            % Make sure daily data are sorted
            [dDaily, ix] = sort(dDaily);
            Xdaily = Xdaily(ix,:);
            
            % ---- Build weekly-averaged daily rates aligned to your weekDates ----
            nWeeks = numel(weekDates);
            nAge   = size(Xdaily,2);
            Xweek  = nan(nWeeks, nAge);
            
            for k = 1:nWeeks
                t0 = weekDates(k);
                if k < nWeeks
                    t1 = weekDates(k+1);     % next window start (exclusive end)
                else
                    % last window: end at last available daily date + 1 day (exclusive)
                    t1 = max(dDaily) + 1;
                end
            
                idx = (dDaily >= t0) & (dDaily < t1);
                if any(idx)
                    Xweek(k,:) = mean(Xdaily(idx,:), 1, 'omitnan');  % daily-averaged rate over that window
                else
                    % no daily data in this window -> leave NaN or fill how you prefer
                    Xweek(k,:) = NaN;
                end
            end
            
            % ---- Output weekly table using SAME weekly dates as calibration ----
            sigmaWeeklyTable = array2table(Xweek, 'VariableNames', sigmaDailyTable.Properties.VariableNames(2:end));
            sigmaWeeklyTable.date = datetime(weekDates, 'ConvertFrom','datenum');
            sigma = movevars(sigmaWeeklyTable, 'date', 'Before', 1);            
        end
  

        % function sigma2 = get.sigma2(obj)
        %     % ---- Inputs ----
        %     sigmaDailyTable = readtable(obj.sigmasHandles(2));  % has sigmaDailyTable.date + age cols
        %     [sigmaDailyTable, isError] = uniqueVaccineDatesFunc(sigmaDailyTable, obj.dataFolder);
        %     if isError
        %         logAndThrowError(obj, 'Daily sigma dates are not unique.');
        %     end
        % 
        %     weekDates = obj.dateNumList(logical(obj.vaccinate));   % datenum weekly list (your calibration list)
        % 
        %     % ---- Convert daily table to arrays ----
        %     dDaily = datenum(sigmaDailyTable.Var1);
        %     Xdaily = table2array(sigmaDailyTable(:,2:end));  % [Ndaily x Nage]
        % 
        %     % Make sure daily data are sorted
        %     [dDaily, ix] = sort(dDaily);
        %     Xdaily = Xdaily(ix,:);
        % 
        %     % ---- Build weekly-averaged daily rates aligned to your weekDates ----
        %     nWeeks = numel(weekDates);
        %     nAge   = size(Xdaily,2);
        %     Xweek  = nan(nWeeks, nAge);
        % 
        %     for k = 1:nWeeks
        %         t0 = weekDates(k);
        %         if k < nWeeks
        %             t1 = weekDates(k+1);     % next window start (exclusive end)
        %         else
        %             % last window: end at last available daily date + 1 day (exclusive)
        %             t1 = max(dDaily) + 1;
        %         end
        % 
        %         idx = (dDaily >= t0) & (dDaily < t1);
        %         if any(idx)
        %             Xweek(k,:) = mean(Xdaily(idx,:), 1, 'omitnan');  % daily-averaged rate over that window
        %         else
        %             % no daily data in this window -> leave NaN or fill how you prefer
        %             Xweek(k,:) = NaN;
        %         end
        %     end
        % 
        %     % ---- Output weekly table using SAME weekly dates as calibration ----
        %     sigmaWeeklyTable = array2table(Xweek, 'VariableNames', sigmaDailyTable.Properties.VariableNames(2:end));
        %     sigmaWeeklyTable.date = datetime(weekDates, 'ConvertFrom','datenum');
        %     sigma2 = movevars(sigmaWeeklyTable, 'date', 'Before', 1);          
        % end        
        % 
        % function sigmaB = get.sigmaB(obj)
        %     % ---- Inputs ----
        %     sigmaDailyTable = readtable(obj.sigmasHandles(3));  % has sigmaDailyTable.date + age cols
        %     [sigmaDailyTable, isError] = uniqueVaccineDatesFunc(sigmaDailyTable, obj.dataFolder);
        %     if isError
        %         logAndThrowError(obj, 'Daily sigma dates are not unique.');
        %     end
        % 
        %     weekDates = obj.dateNumList(logical(obj.vaccinate));   % datenum weekly list (your calibration list)
        % 
        %     % ---- Convert daily table to arrays ----
        %     dDaily = datenum(sigmaDailyTable.Var1);
        %     Xdaily = table2array(sigmaDailyTable(:,2:end));  % [Ndaily x Nage]
        % 
        %     % Make sure daily data are sorted
        %     [dDaily, ix] = sort(dDaily);
        %     Xdaily = Xdaily(ix,:);
        % 
        %     % ---- Build weekly-averaged daily rates aligned to your weekDates ----
        %     nWeeks = numel(weekDates);
        %     nAge   = size(Xdaily,2);
        %     Xweek  = nan(nWeeks, nAge);
        % 
        %     for k = 1:nWeeks
        %         t0 = weekDates(k);
        %         if k < nWeeks
        %             t1 = weekDates(k+1);     % next window start (exclusive end)
        %         else
        %             % last window: end at last available daily date + 1 day (exclusive)
        %             t1 = max(dDaily) + 1;
        %         end
        % 
        %         idx = (dDaily >= t0) & (dDaily < t1);
        %         if any(idx)
        %             Xweek(k,:) = mean(Xdaily(idx,:), 1, 'omitnan');  % daily-averaged rate over that window
        %         else
        %             % no daily data in this window -> leave NaN or fill how you prefer
        %             Xweek(k,:) = NaN;
        %         end
        %     end
        % 
        %     % ---- Output weekly table using SAME weekly dates as calibration ----
        %     sigmaWeeklyTable = array2table(Xweek, 'VariableNames', sigmaDailyTable.Properties.VariableNames(2:end));
        %     sigmaWeeklyTable.date = datetime(weekDates, 'ConvertFrom','datenum');
        %     sigmaB = movevars(sigmaWeeklyTable, 'date', 'Before', 1);          
        % end           

         function vaccineCoverage = get.vaccineCoverage(obj)
             % Read the vaccine coverage table from a file handle
             vaccineCoverageTable = readtable(obj.vaccineCoverageHandle);

             % Ensure that dates in the vaccine coverage table are unique
             [vaccineCoverageTable, isError] = uniqueVaccineDatesFunc(vaccineCoverageTable, obj.dataFolder);
             
             % Log and throw an error if dates are not unique as expected
             if isError
                 logAndThrowError(obj, 'Vaccine coverage dates are not unique.');
             end

% %              % Check that the dates are contiguous (no gaps)
% %              isError = vaccineCoverageNotContiguousFunc(vaccineCoverageTable);
% %              % Log and throw an error if dates are not contiguous
% %              if isError
% %                  logAndThrowError(obj, 'Vaccine coverage dates are not contiguous.');
% %              end

             % Filter for specific dates based on vaccination data availability
             vaccineCoverageMatch = vaccineCoverageTable(ismember(...
                 datenum(vaccineCoverageTable.report_date),...
                 obj.dateNumList(logical(obj.vaccinate))), :);

             if obj.vaccineCoverageLabel==0

                 % Calculate vaccine coverage as a percentage of the total population
                 vaccineCoverage = [...
                     vaccineCoverageMatch.total_individuals_at_least_one / obj.totalPopulation, ...
                     vaccineCoverageMatch.total_individuals_fully_vaccinated / obj.totalPopulation, ...
                     vaccineCoverageMatch.total_individuals_3doses / obj.totalPopulation] * 100;
    
                 % Replace NaN values with zero in vaccine coverage data
                 vaccineCoverage(isnan(vaccineCoverage)) = 0;
             else
                 vaccineCoverage = [...
                     vaccineCoverageMatch.coverage_at_least_one , ...
                     vaccineCoverageMatch.coverage_fully_vaccinated, ...
                     vaccineCoverageMatch.coverage_3doses ];
             end

             % Apply a decay function to the booster data if the flag is set
             if obj.boosterDecayFlag
                 vaccineCoverage(:,:) = calculateBoosterDecay(obj, vaccineCoverage, 0.1);
             end

             % Calculate the change in coverage for the booster data and prepend a zero
             obj.deltaCoverageBoost = diff(vaccineCoverage(:, 3));
             obj.deltaCoverageBoost = [0; obj.deltaCoverageBoost];
         end

        function tmax = get.tmax(se)
            tmax = (se.dateNumList(end)-se.dateNumList(1));
        end
        function N = get.N(obj)
            N = length(obj.age);
        end
        function UStotal = get.UStotal(obj)
            UStotal = obj.USmales+obj.USfemales;
        end
        function PerFemales = get.PerFemales(obj)
            PerFemales = obj.USfemales./obj.UStotal;
        end
        %         function f = get.f(obj)
        %             f = obj.fRaw.*obj.PerFemales;
        %         end
        function gamma = get.gamma(obj)
            gamma(1,:) = NaN(1,obj.N);
            gamma(2,:) = obj.gamma2*ones(1,obj.N);
            gamma(3,:) = obj.gamma3*ones(1,obj.N);
            gamma(4,:) = obj.gamma4*ones(1,obj.N);
        end
        function beta = get.beta(obj)
            beta(1,:) = zeros(1,obj.N);
            beta(2,:) = 0.5.*obj.beta3*ones(1,obj.N);
            beta(3,:) = obj.beta3*ones(1,obj.N);
            beta(4,:) = 0.1.*obj.beta3*ones(1,obj.N);
        end
        function p2 = get.p2(obj) % % infection fraction from S2 to I2,I3,I4
            p2 = obj.p1;
        end
        function p3 = get.p3(obj) % infection fraction from S3 to I2,I3,I4
            p3 = obj.p1;
            p3(1,:) = 0;
        end
        function delta = get.delta(obj)
            delta = ones(4,obj.N);
            delta(1:2,1:obj.N) = 0;
            delta(3,1:obj.N) = obj.delta3;
            delta(4,1:obj.N) = obj.delta4;
        end
        function kappa = get.kappa(obj)
            kappa2 = 1/obj.kappaZero;
            kappa3 = 1/obj.kappaZero;
            kappa4 = 1/obj.kappaZero;
            kappa = ones(4,obj.N);
            kappa(1,:) = obj.kappa1*kappa(1,:);
            kappa(2,:) = kappa2*kappa(2,:);
            kappa(3,:) = kappa3*kappa(3,:);
            kappa(4,:) = kappa4*kappa(4,:);
        end
        function tau = get.tau(obj)
            tau = NaN(1,obj.N);
            for i = 1:obj.N
                tau(i) = (obj.mu(i)+obj.q)/...
                    (exp((obj.mu(i)+obj.q)*365*obj.age_interval(i))-1);
                tau(i) = 1/(5*365); % fixed rate overriding the above formula
            end
            tau(end) = 0;
        end
        function alpha = get.alpha(obj)
            alpha = obj.alphaInput'.*ones(1,obj.N);
        end
        function pastPredictionDate = get.pastPredictionDate(obj)
            pastPredictionDate = obj.hospDates(end);
        end

        %% Data plot
        function plotData(obj)
            figure
            subplot(3,1,1)
            plot(obj.incidenceDates,obj.incidenceCases,'+','MarkerSize',1)
            hold on
            plot(obj.hospDates,obj.hospCases,'o','MarkerSize',1)
            datetick('x','mmm-yy','keepticks')
            ylabel("Daily Cases")
            vocLines(obj)
            legend("Incidence","Hospitalizations")

            subplot(3,1,2)
            plot(obj.incidenceDates,obj.incidenceCases,'+','MarkerSize',1)
            hold on
            plot(obj.hospDates,obj.hospCases,'o','MarkerSize',1)
            datetick('x','mmm-yy','keepticks')
            ylabel("Daily Cases")
            ylim([0,3000])
            vocLines(obj)
            legend("Incidence","Hospitalizations") 

            % Find common dates and their indices between the two datasets
            [Dates_c,ind_incidence_c]=intersect(obj.incidenceDates,obj.hospDates);
            [~,ind_hosp_c]=intersect(obj.hospDates,obj.incidenceDates);

            subplot(3,1,3)
            plot(Dates_c,obj.hospCases(ind_hosp_c)./obj.incidenceCases(ind_incidence_c),'+','MarkerSize',1,'Color',[0.4940 0.1840 0.5560])
            datetick('x','mmm-yy','keepticks')
            ylabel("Fraction Hospitalized")
            vocLines(obj)
        end

        
        %% parameter plots
        function boosterCoveragePlot(obj, ~)
            if nargin < 2
                    figure;
                    ax = gca;
            end
            boosterDosePlot(obj);
            vertLines(obj);
            vocLines(obj);

            saveBool = true;
            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', 'booster dose');
                print(fullFileName,'-dpng','-r300')
            end
        end
        %% simulation
        %% three methods to run simulation base simulation

        function simulate(obj,simulateWindow) %run simulation (aka create yout)
            arguments
                obj
                simulateWindow = [1 length(obj.dateList)];
            end
            %% not working short circit
            %simulateWindow = [1 simulateWindow(2)];
            %%
            if obj.cut_Omicron==1
                obj.contactModifers=obj.contactModifers(1:length(obj.dateList),:);
            end
            if length(obj.contactModifers)~=length(obj.dateNumList)
                disp('Contact Modifier length error: Are the prediction values set?')
            end
            mw = timeBreaksFunc(obj.dateList); % timeBreaks (difference in days from day0)
            [~,ps]=find(mw<=simulateWindow(1),1,"last"); % start timeBreak index
            [~,pe]=find(mw<=simulateWindow(2),1,"last"); % end timeBreak index
            ps = max(ps,3); %avoid negative first window
            %startIndex = (obj.dateNumList(mw(ps)+1)-obj.dateNumList(1))+1;
            if simulateWindow(1)==1%||mw(ps)<mw(2)
                if simulateWindow(2)<length(obj.dateList)
                    simulateWindow = [1,mw(pe)-1];
                end
                startIndex = 1;
                yint = immunoEpiBoostYINT(obj, obj.population);
            else
                startIndex = (obj.dateNumList(mw(ps-1))-obj.dateNumList(1))+1;
                yint = obj.yout(startIndex,:)';%time starts at 0
                simulateWindow = [mw(ps-1),mw(pe)-1];
            end
            [timeOut,yOutWindow,~]=...
                RunImmunoEpiFast_Boost(obj,...
                obj.kval, ...
                simulateWindow,...
                yint);
            obj.yout = [obj.yout(1:startIndex-1,:);yOutWindow];%[yOutStart;yOutWindow];
            obj.time = [(0:startIndex-2)';timeOut+startIndex-1];%+timestart(end);
        end

        %run simulation and return [time, yout] 
        function [time,yout] = simulateReturn(obj,simulateWindow) %run simulation (aka create yout)
            arguments
                obj
                simulateWindow = [1 length(obj.dateList)];
            end
            %% not working short circit
            %simulateWindow = [1 simulateWindow(2)];
            %%
            if length(obj.contactModifers)~=length(obj.dateNumList)
                disp('Contact Modifier length error: Are the prediction values set?')
            end
            mw = timeBreaksFunc(obj.dateList);
            [~,ps]=find(mw<=simulateWindow(1),1,"last");
            [~,pe]=find(mw>=simulateWindow(2),1,"first");
            ps = max(ps,3);%avoid negative first window
            %startIndex = (obj.dateNumList(mw(ps)+1)-obj.dateNumList(1))+1;
            if simulateWindow(1)==1%||mw(ps)<mw(2)
                if simulateWindow(2)<length(obj.dateList)
                    simulateWindow = [1,mw(pe)-1];
                end
                startIndex = 1;
                yint = immunoEpiBoostYINT(obj, obj.population);
            else
                startIndex = (obj.dateNumList(mw(ps-1))-obj.dateNumList(1))+1;
                yint = obj.yout(startIndex,:)';%time starts at 0
                simulateWindow = [mw(ps-1),mw(pe)-1];
            end
            [timeOut,yOutWindow,~]=...
                RunImmunoEpiFast_Boost(obj,...
                obj.kval, ...
                simulateWindow,...
                yint);
            yout = [obj.yout(1:startIndex-1,:);yOutWindow];%[yOutStart;yOutWindow];
            time = [(0:startIndex-2)';timeOut+startIndex-1];%+timestart(end);
        end

        % run simulation with window speed enhancement for fiting
        function [time,yout] = fitSimulate(obj,simulateWindow) %run simulation (aka create yout)
            arguments
                obj
                simulateWindow = [1 length(obj.dateList)];
            end
            %obj.simulate
            % not working short circit
            % simulateWindow = [1 simulateWindow(2)];
            %%
            %             if length(obj.contactModifers)~=length(obj.dateNumList)
            %                 disp('Contact Modifier length error: Are the prediction values set?')
            %             end
            %             mw = timeBreaksFunc(obj.dateList);
            %              [~,ps]=find(mw<=simulateWindow(1),1,"last");
            %              [~,pe]=find(mw>=simulateWindow(2),1,"first");
            %            ps = max(ps,3);%avoid negative first window
            %startIndex = (obj.dateNumList(mw(ps)+1)-obj.dateNumList(1))+1;
            if simulateWindow(1)==1%||mw(ps)<mw(2)
                %                 if simulateWindow(2)<length(obj.dateList)
                %                 simulateWindow = [1,mw(pe)-1];
                %                 end
                startIndex = 1;
                yint = immunoEpiBoostYINT(obj, obj.population);
            else
                startIndex = (obj.dateNumList(simulateWindow(1))-obj.dateNumList(1))+1;
                if startIndex>size(obj.yout,1)
                    disp("Error")
                end
                yint = obj.yout(startIndex,:)';%time starts at 0
                %simulateWindow = [mw(ps-1),mw(pe)-1];
            end
            [timeOut,yOutWindow,~]=...
                RunImmunoEpiFast_Boost(obj,...
                obj.kval, ...
                simulateWindow,...
                yint);
            yout = [obj.yout(1:startIndex-1,:);yOutWindow];%[yOutStart;yOutWindow];
            time = [(0:startIndex-2)';timeOut+startIndex-1];%+timestart(end);
        end

        % run simulation with seroprevalence details
                function simulateDetail(obj,simulateWindow) %run simulation (aka create yout)
            arguments
                obj
                simulateWindow = [1 length(obj.dateList)];
            end
            %% not working short circit
            %simulateWindow = [1 simulateWindow(2)];
            %%
            if length(obj.contactModifers)~=length(obj.dateNumList)
                disp('Contact Modifier length error: Are the prediction values set?')
            end
            mw = timeBreaksFunc(obj.dateList);
            [~,ps]=find(mw<=simulateWindow(1),1,"last");
            [~,pe]=find(mw>=simulateWindow(2),1,"first");
            ps = max(ps,3);%avoid negative first window
            %startIndex = (obj.dateNumList(mw(ps)+1)-obj.dateNumList(1))+1;
            if simulateWindow(1)==1%||mw(ps)<mw(2)
                if simulateWindow(2)<length(obj.dateList)
                    simulateWindow = [1,mw(pe)-1];
                end
                startIndex = 1;
                yint = immunoEpiBoostYINTDetail(obj, obj.population);
            else
                startIndex = (obj.dateNumList(mw(ps-1))-obj.dateNumList(1))+1;
                yint = obj.yout(startIndex,:)';%time starts at 0
                simulateWindow = [mw(ps-1),mw(pe)-1];
            end
            [timeOut,yOutWindow,~]=...
                RunImmunoEpiFast_Boost_Detail(obj,...
                obj.kval, ...
                simulateWindow,...
                yint);
            obj.youtDetail = [obj.yout(1:startIndex-1,:);yOutWindow];%[yOutStart;yOutWindow];
            obj.timeDetail = [(0:startIndex-2)';timeOut+startIndex-1];%+timestart(end);
        end

        %% utilities

        function loadKvalue(obj) %reload kvalue from text file
            kvalTable = readtable(obj.kValueHandle);
            obj.kval = kvalTable.Var1;

            lenKval = length(obj.kval);
            lenDateList = length(obj.dateList);

            if lenKval ~= lenDateList
                disp("Error: kValue List not equal length of date list");
                disp(['Length of kval: ', num2str(lenKval)]);
                disp(['Length of dateList: ', num2str(lenDateList)]);
            end

            % Write to simulation Log file

            logFile = fopen(obj.errorLogFilePath, 'a');
            if logFile == -1
                error('Unable to open the log file.');
            end
            % Write the debug information to the file
            fprintf(logFile, '%s - obj.loadKvalue\n', datetime('now'));
            fprintf(logFile, '- Error: kValue List not equal length of date list\n');
            fprintf(logFile, '- Length of kval: %s , num2str(lenKval)\n', num2str(lenKval));
            fprintf(logFile, '- Length of dateList: %s\n', num2str(lenDateList));

            % Close the log file
            fclose(logFile);

            disp(obj.kValueHandle);
            obj.simulate;
        end


        function r = post(obj) %transform results into human readable form
            r = postProcess(obj.yout,16);
        end

        function window = mitigationIndex2Window(obj,index)
            mitagationWindows = pickMitagationWindows(obj,0,0);
            first = find(mitagationWindows>=index,1,"first");
            last = find(mitagationWindows<=index,1,"last");
            window = first;
            if first~=last
                disp(['range: ',num2str(first),':',num2str(last)])
            end
        end

        function kvalueTail(obj,type)
            arguments
                obj
                type = 'mean'
            end
            meanLength = 7*6;

            if obj.hospital ==1
                listSplitIndex = find(obj.dateNumList<=obj.hospDates(end),1,"last");
            else
                listSplitIndex = find(obj.dateNumList<=obj.incidenceDates(end),1,"last");
            end

            meanLength = min(listSplitIndex,meanLength);

            if strcmp(type,'mean')
                kvalueTailMeanValue = mean(obj.kval(listSplitIndex-meanLength:listSplitIndex));
                kvalueTail = kvalueTailMeanValue.*ones(length(obj.kval(listSplitIndex+1:end)),1);
            elseif strcmp(type,'lastValue')
                kvalueTailValue = obj.kval(listSplitIndex);
                kvalueTail = kvalueTailValue.*ones(length(obj.kval(listSplitIndex+1:end)),1);
            else
                disp('no change')
                kvalueTail = obj.kval(listSplitIndex+1:end);
            end
            obj.kval = [obj.kval(1:listSplitIndex);kvalueTail];
        end

        %% prediction functions
        function predictionKval(obj,opt)
            predictedKval = predictionKvalFunc(obj,opt);

            %             obj.kval = predictedKval(:,2);
            %             [~,obj.youtUB] = obj.simulateReturn;
            %
            %             obj.kval = predictedKval(:,3);
            %             [~,obj.youtLB] = obj.simulateReturn;

            obj.kval = predictedKval(:,1);
            obj.simulate
        end

        function predictionOmega(obj,opt,waveFactor)
            arguments
                obj
                opt
                waveFactor = 12;
            end

            obj.omega = predictionOmegaFunc(obj,opt,waveFactor);
            obj.simulate
        end

        function predictionContacts(obj,opt)
            predictionIndex = obj.hospDates(end)<=obj.dateNumList;
            contactModifierOther = ones(sum(predictionIndex),1);
            if strcmp(opt,'10')
                obj.contactModifers(predictionIndex,2) = 5.*contactModifierOther;
            elseif strcmp(opt,'20')
                obj.contactModifers(predictionIndex,2) = 6.*contactModifierOther;
            elseif strcmp(opt,'30')
                obj.contactModifers(predictionIndex,2) = 3.*contactModifierOther;
            else
                disp('no change')
            end
        end

        %% fit
        function fit(obj)
            % if obj.hospital==2
            %     asiimFitIII(obj,obj)
            % else 
            %     asiimFitII(obj,obj)
            % end

            %asiimFitII(obj)
            asiimFitIII(obj)

            %obj.kvalueTail
            %writematrix(obj.kval,obj.kValueHandle);
        end

        function ret = get.fitCurvature(obj)
            dataNS = obj.hospCases;
            dateNS = obj.hospDates;
            fitProp.fitProportions = obj.fitProportions;
            endWindow = find(obj.dateNumList<obj.hospDates(end),1,'last');
            objective= forceFitObjFastFunc(obj.kval,obj,dataNS,dateNS,...
                1,endWindow,...
                true,1, 1, fitProp);
            curvaturePenaltyCalc = objective*obj.curvaturePenaltyPercent;
            ret = sum(curvaturePenaltyCalc*diff(obj.kval,2).^2)/objective;
        end
        function error = get.error(obj)
            dataNS = obj.hospCases;
            dateNS = obj.hospDates;
            fitProp.fitProportions = obj.fitProportions;
            endWindow = find(obj.dateNumList<obj.hospDates(end),1,'last');
            error = forceFitObjErrorOutFunc(obj.kval,obj,dataNS,dateNS,...
                1,endWindow,...
                true,1, 1, fitProp);
        end

        %% fit diagonostics
        function mse = get.mse(obj)
            mse = obj.error.^2;
        end
        function percentError(obj,yLimits)
            arguments
                obj
                yLimits = 'auto';%[3e-3,10e-3]
            end
            percentError = smooth(obj.error,'lowess')./smooth(obj.hospCases,'lowess');
            figure()
            hold on
            bar(obj.hospDates,percentError);
            %             yyaxis left
            %             bar(obj.fitError);
            hold off
            datetick('x','mmm-yy','keepticks')

            h=gcf;
            set(h,'PaperType','uslegal');
            set(h,'PaperPositionMode','auto');
            set(h,'PaperOrientation','landscape');
            set(h,'Position',pp.p);
            ylabel('$\%$ error','interpreter',...
                'latex','FontSize',14)
            xlabel('date','interpreter',...
                'latex','FontSize',14);
            xlim('padded')
            ylim(yLimits)


            fullFileName = fullfile(obj.folder,'/results/', 'percent error');
            print(fullfile(obj.experimentDir(),'\',fullFileName),'-dpng','-r300')
        end

        function kvalVar = get.kvalVar(obj)
            kvalVar = kvalVarFunc(obj);
        end

        function AIC = AIC(obj)
            AIC = aicFunc(sum(obj.forceFit),length(obj.hospCases),obj.mse);
        end
        function BIC = BIC(obj)
            BIC = bicFunc(sum(obj.forceFit),length(obj.hospCases),obj.mse);
        end
        function normalTest = get.normalTest(obj)
            %Kolmogorov-Smirnov for unknown mean and var
            [h,~,normalTest,critical]=lillietest(obj.error,'MCTol',1e-4);
            disp(['reject',num2str(h),' kstat:', num2str(normalTest),' critical: ',num2str(critical)])
        end

        %% results

        function I4 = I4(obj)
            rp = postProcess(obj.yout);
            I4 = sum(rp(4).I+rp(8).I,2);
        end
        function ret = incidence(obj)
            rp = postProcess(obj.yout);
            ret = convertDaily(rp);
        end

        function saveActive(obj, name)
            %%saveActive(obj, name) Save the active figure
            arguments
                obj
                name = 'custom figure';
            end
            fullFileName = fullfile(obj.folder,'/results/', name);
            print(fullfile(obj.experimentDir(),'\',fullFileName),'-dpng','-r300')
        end

        function fitPlot(obj,fast,fitWindows,screen,title, savePlot, subPlot)
            %% fit plot
            arguments
                obj
                fast = true;
                fitWindows = [0,0];
                screen = false;
                title = 'detailed fit';
                savePlot = true;
                subPlot = 'no';
            end
            if strcmp(subPlot,'no')
                figure;
                ax = gca;
            end

            pp = plotProp;
            if screen
                pp.plotSize = 500;
                pp.fontSize = 12;
            end
            
            
            %ret = 
            detailedFitPlot(obj,pp,fitWindows,fast);

            if savePlot
            fullFileName = fullfile(obj.folder, title);
            print(fullFileName,'-dpng','-r300')
            end
        end

        function sensitivityFitPlot(obj,pp,fitWindows,fast,title)
            %% fit plot
            arguments
                obj
                pp = plotProp;
                fitWindows = [0,0];
                fast = true;
                title = 'detailed fit';
            end
            screen = true;
            %pp = plotProp;
            if screen
                pp.plotSize = 500;
                pp.fontSize = 12;
            end
            sensitivityPlot(obj,pp,fitWindows,fast)
            fullFileName = fullfile(obj.folder,'/results/', title);
            print(fullfile(obj.experimentDir(),'\',fullFileName),'-dpng','-r300')
        end
        function contactReduction(obj,saveBool,title,visibilityFlag)
            arguments
                obj
                saveBool = true;
                title = 'contact reduction';
                visibilityFlag = 'on';
            end
            
            contactPlot(obj,visibilityFlag)

            if saveBool
                fullFileName = fullfile(obj.folder, title);
                print(fullfile(obj.experimentDir,fullFileName),'-dpng','-r300')
            end
        end

        function kValuePlot(obj, saveBool)
            arguments
                obj
                saveBool = true;
            end
            kValueBar(obj)

            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', 'k value');
                print(fullFileName,'-dpng','-r300')
            end
        end

        function seroprevalence(obj,saveBool,subplotBool,legendBool)
            arguments
                obj
                saveBool = true;
                subplotBool = false;
                legendBool = true;
            end
            pp = plotProp; %create a plotProperty object from class
            %pp.fileName = ['Seroprevalence','-',datestr(obj.stopDate)];
            pp.showLegend = legendBool;
            pp.subPlot = subplotBool;

            if obj.hospital ==1
                data.dates = obj.hospDates;
            else
                data.dates = obj.incidenceDates;
            end

            %figure()
            %resutsPlotsSeroprevalence(obj,obj,data,pp)
            resutsPlotsSeroprevalence2(obj,4:16)

            %vertLines(obj)

            %title(obj.region,'interpreter','latex','FontSize',pp.fontSize+2)
            ylabel('Population Proportion','interpreter',...
                'latex','FontSize',pp.fontSize-2)
            %xlabel('time','interpreter','latex','FontSize',pp.fontSize+2)

            %h=gcf;
            %set(h,'PaperType','uslegal');
            %set(h,'PaperPositionMode','auto');
            %set(h,'PaperOrientation','landscape');
            %set(h,'Position',pp.position);

            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', 'seroprevalence');
                print(fullFileName,'-dpng','-r300')
            end
        end



        %%
        function immunity(obj,saveBool,title)
            arguments
                obj
                saveBool = true;
                title='';
            end
            pp = plotProp;

            pp.percent = true;

            %save flag for subfunction (not at end of this function!)
            pp.save = false;
            pp.showLegend = true;

            pp.fileName = ['','SeroprevalenceImunity',...
                '-',datestr(obj.stopDate),'-',title];
            pp.region = '';

            if obj.hospital ==1
                data.dates = obj.hospDates;
            else
                data.dates = obj.incidenceDates;
            end

            resutsPlotsSeroprevalenceImmunityBoost(obj,obj,data,pp)
            %title(obj.shortTitle,'interpreter','latex','FontSize',pp.fontSize+2)
            ylabel('Population Proportion','interpreter',...
                'latex','FontSize',pp.fontSize-2)
            %xlabel('time','interpreter','latex','FontSize',pp.fontSize+2)

            % h=gcf;
            % set(h,'PaperType','uslegal');
            % set(h,'PaperPositionMode','auto');
            % set(h,'PaperOrientation','landscape');
            % %set(h,'Position',[-1800 70 1600 900]);
            % set(h,'Position',pp.position);

            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', 'immunity');
                print(fullFileName,'-dpng','-r300')
            end
        end

        function seroFitPlot(obj,saveBool)
            arguments
                obj
                saveBool = true;
            end

            pp = plotProp;

            pp.ageRange = 1:16;
            pp.showTitle = 0;
            pp.save = 0;
            pp.visible = 'on';
            pp.subPlot = 1;
            pp.fontSize = 10;
            pp.percent = 1;
            pp.fontSize = 22;

            pp.fileName = ['SeroprevalenceNatural',...
                '-',datestr(obj.stopDate)];
            pp.region = '';

            data.dates = obj.hospDates;

            figure()
            resutsPlotsSeroprevalenceNatural(obj,obj,data,pp)
            title(obj.shortTitle,'interpreter','latex','FontSize',pp.fontSize+2)
            ylabel('population proportion','interpreter',...
                'latex','FontSize',pp.fontSize+2)
            xlabel('time','interpreter','latex','FontSize',pp.fontSize+2)

            obj.vocLines
            obj.vertLines

            h=gcf;
            set(h,'PaperType','uslegal');
            set(h,'PaperPositionMode','auto');
            set(h,'PaperOrientation','landscape');
            set(h,'Position',[-1800 70 1600 900]);

            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', 'seroFit');
                print(fullfile(obj.experimentDir(),'\',fullFileName),'-dpng','-r300')
            end
        end

        function immuneBurden(obj, I4flag, saveBool)
            arguments
                obj
                I4flag = true;
                saveBool = true;
            end
            pp = plotProp;
            if I4flag
                pp.fileName = 'immuneByAgeBurden-I4';
                pp.I4Only = true;
            else
                pp.fileName = 'immuneByAgeBurden';
                pp.I4Only = false;
            end

            pp.plotSize = 400;

            figure()
            resultsImmuneBurdenByAgeSubPlot(obj,obj,pp);
            h=gcf;
            set(h,'PaperType','uslegal');
            set(h,'PaperPositionMode','auto');
            set(h,'PaperOrientation','landscape');
            set(h,'Position',pp.position);

            if saveBool
                fullFileName = fullfile(obj.folder,'/results/', pp.fileName);
                print(fullfile(obj.experimentDir(),'\',fullFileName),'-dpng','-r300')
            end
        end

        %% validation methods
        function compartments(obj,ageGroups)
            arguments
                obj
                ageGroups = 1:16;
            end
            if length(obj.yout)<2
                disp('error: re-run simulation!')
            else
                %compartmentsPlot(obj,ageGroups);
                %compartmentsPlot_v2(obj, ageGroups);
                compartmentsPlot_v3(obj, ageGroups);
                compartmentsPlot_Vaccinations(obj, ageGroups);
            end
        end

        %% Diagnostic Plots
        function plotSimulationSetup(obj)
            plotSimulationSetupFunc(obj)
        end

        function plotPax(obj, ax)
                if nargin < 2
                    figure;
                    ax = gca;
                end
            plotPaxVsDate(obj,ax)
            vertLines(obj)
            vocLines(obj)
        end

        function plotOmega(obj, ax)
               if nargin < 2
                    figure;
                    ax = gca;
                end
            plotOmegaVsDate(obj,ax)
            vertLines(obj)
            vocLines(obj)
        end


        %% Overplots
        function I4OverPlot(obj)
            rp = postProcess(obj.yout);
            I4 = sum(rp(4).I+rp(8).I,2)/sum(obj.population);
            pp = plotProp;

            hold on
            pl1 = plot(obj.time+obj.dateNumList(1),1-I4,'k');
            pl2 = plot(obj.time+obj.dateNumList(1),1-I4,'r');
            pl2.LineStyle = '-.';
            pl1.LineWidth = pp.lineWidth+.3;
            pl2.LineWidth = pp.lineWidth+.1;
            hold off
        end

        function overPlot(obj)
            %overlay kval, and vaccine coverage
            %plot range [0,1]
            hold on
            overPlots(obj)
            hold off
        end

        function vaccineOverPlot(obj)
            vaccineOverPlotFunc(obj)
        end

        function vertLines(obj)
            hold on
            if obj.hospital ==1
                verticalLines(obj,obj.hospDates(end))
            else
                verticalLines(obj,obj.incidenceDates(end))
            end
            hold off
        end

        function vocLines(obj)
            hold on
            xl(1) = xline(obj.alphaWaveStart,'-.');
            xl(1).Label = 'Alpha';
            xl(2) = xline(obj.deltaWaveStart,'-.');
            xl(2).Label = 'Delta';
            xl(3) = xline(obj.omicronWaveStart,'-.');
            xl(3).Label = 'Omicron';
            hold off
        end

        function citfOverPlot(obj)
            rawData = readtable([obj.experimentDir(),'\data\citf-seroprevalence.xlsx']);

            indexProvince = strcmpi(rawData.geo,'ON');
            indexN = strcmpi(rawData.ab_estimate,'N');
            indexS = strcmpi(rawData.ab_estimate,'S');

            pctNMean = rawData.pct_mean(logical(indexProvince.*indexN));
            datesN =  datenum(rawData.endofweek(logical(indexProvince.*indexN)));

            pctSMean = rawData.pct_mean(logical(indexProvince.*indexS));
            datesS =  datenum(rawData.endofweek(logical(indexProvince.*indexS)));

            hold on
            p(1) = plot(datesN,pctNMean);
            p(4) = plot(datesS(2:end),pctSMean(2:end));
            hold off
            %             lh = legend({'CITF Seroprevalence Estimates: anti N', ...
            %                 'CITF Seroprevalence Estimates: anti S'}, ...
            %                 'Location','northwest','Orientation','vertical', ...
            %                 'FontSize',12)

            for i = [1,4]
                p(i).LineWidth = 1.4;
            end
            p(1).LineStyle = '-.';
            p(1).Color = 'r';
            p(4).LineStyle = '-.';
            p(4).Color = 'b';
        end

        %% prediction
function modifyKvalue(obj, percentage, numWeeks)
% Modifies kValue based on a percentage of the average
% kValue over a given number of weeks prior to the last
% hospitalization date.
% 
%   obj        - Object containing date and kValue data.
%   percentage - The percentage by which to adjust future kValues.
%                Default is 0.
%   numWeeks   - The number of weeks used to calculate the average kValue.
%                Default is 6 weeks.
%
% Example usage: modifyKvalue(obj, 10, 6)

                arguments
                obj
                percentage = 0;
                numWeeks = 6;
                end

    % Find the index of the last date in obj.hospDates
    lastHospDateIdx = find(obj.dateNumList == obj.hospDates(end), 1);

    if isempty(lastHospDateIdx)
        disp('Error: Last hospDate not found in dateList');
        return;
    end

    daysPerWeek = 7;

    % Calculate the start index for numWeeks weeks prior
    startIdx = max(1, lastHospDateIdx - numWeeks * daysPerWeek);

    % Extract the kValue subset for numWeeks weeks prior
    sixWeeksKvals = obj.kval(startIdx:lastHospDateIdx);

    % Calculate the average kValue for the six-week period
    avgKval = mean(sixWeeksKvals);

    % Modify kValue for the remaining dates after last hospDate
    obj.kval(lastHospDateIdx+1:end) = obj.kval(lastHospDateIdx+1:end) + avgKval * (percentage / 100);
end

%% Calculate vaccination rates over time and plot them

function calculateSigmas(obj)
    arguments
    obj
    end

    [numDates,numDose]=size(obj.vaccineCoverage);
    vacStep = 0; %Tracks how far the vaccination schedule has progressed.
    
    sigma = zeros(4,16,numDose,length(obj.dateNumList));
    p=setParamters(obj,sigma);

    for index=1:length(obj.dateNumList)
        maxCoverage=obj.vaccineCoverage;

        if obj.vaccinate(index) && vacStep<numDates-1
            vacStep = vacStep + 1;
            %caculate vaccination rate sigma
            [sigma(:,:,1,index),~] = sigmaCalculationAgeB(...
                vacStep,...
                p,...
                obj.yout,...%youtSave,...
                obj.newVacTarget,...
                obj.mitagationWindowLengthShort,...
                maxCoverage, ...
                obj.population);
            
            %caculate vaccination rate sigma2
            [sigma(:,:,2,index),~] = sigmaCalculationTwoB(...
                vacStep,...
                p,...
                obj.yout,...%youtSave,...
                obj.mitagationWindowLengthShort,...
                maxCoverage(:,2), ...
                obj.population);
              
            [sigma(:,:,3,index),~] = sigmaCalculationBoostB(...
                vacStep,...
                p,...
                obj.yout,...%youtSave,...
                obj.newVacTarget,...
                maxCoverage(:,3),...
                obj.mitagationWindowLengthShort,...
                obj.deltaCoverageBoost, ...%delta change
                obj.population);
    
             if index > 200
                 %bug hunt
                 %disp('bug hunt active! stuff set to 0')
                 %disp(['index: ',num2str(index)])
                 %disp(p.sigma)
                 p.sigma = p.sigma*.0;
                 p.sigma2 = p.sigma2*.0;
                 p.sigmaB = p.sigmaB;
             end
        end

    end

end

    end
end

