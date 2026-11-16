# Italy COVID-19 Application

This directory contains the MATLAB code for applying the **Age-Structured Immunity and Infection Model (ASIIM)** to COVID-19 in Italy. The application uses the common functions in [`../functions/`](../functions/) and Italy-specific epidemiological data and model settings.

The example fits the model to **reported COVID-19 incidence and hospital admissions** using a rolling-window genetic algorithm to estimate a time-varying transmission factor, `kappa(t)`.

## Directory structure

```text
ASIIM/
├── functions/                             Shared MATLAB functions
├── Italy/
│   ├── example_Italy.m                    Runs one calibration and plots results
│   ├── ItalySetup.m    	           Loads Italy data and model settings
│   ├── data/                              Required Italy input files (add separately)
│   ├── results/                           Generated calibration outputs
│   └── README.md
└── README.md
```

## Requirements

- MATLAB (the original COVID-ASIIM package specified R2021a or later; compatibility has not been verified across releases)
- **Global Optimization Toolbox** for genetic algorithm calibration (`ga`)
- **Statistics and Machine Learning Toolbox** for functions such as `normrnd`
- **Parallel Computing Toolbox** if `sim.useParallel = true`; change this to `false` to disable parallel optimization

## Run the Italy example

1. Clone or download the [ASIIM repository](https://github.com/MarwaTuffaha/ASIIM).
2. Place the required Italy input files inside `Italy/data/`.
3. Set the MATLAB **Current Folder** to `ASIIM/Italy`.
4. Run:

```matlab
example_Italy
```

The example adds `../functions/` to the MATLAB path and creates a simulation object:

```matlab
sim = ItalySetup();
sim.fit
sim.simulate
sim.fitPlot
```

**Each run performs one calibration**, not the 20-calibration averaging block present in the older `example_Italy.m` script. The example uses the original fitting parameters from that script, including fitting windows 50 through 102. The setup file contains the Italy-specific parameter values and input handling.

## Input data

The `data/` directory must contain the following files referenced by `ItalySetup.m`:

| Category | Input files |
| --- | --- |
| Epidemiological observations | `dateList.txt`, `incidenceData.txt`, `OWD_hospAdmissionData.csv` |
| Population and immunity | `population2020.txt`, `p1_Italy.txt`, `omega.txt` |
| Initial transmission factor | `kvalue.csv` |
| Contact modifiers | `contactModifiers_weekly.txt` |
| Contact matrices | `italy_contactHome.txt`, `italy_contactWork.txt`, `italy_contactOtherPlaces.txt`, `italy_contactSchool.txt` |
| Contact perturbation matrices | `contactSF1.txt`–`contactSF3.txt`, `contactWF1.txt`–`contactWF3.txt`, `contactOF1.txt`–`contactOF3.txt` |
| Incidence and hospitalization fitting fractions | `fitPropInc_Var_July2020.csv`, `fitPropHosp_Var_July2020.csv` |
| Vaccination | `vaccinate.txt`, `vaccineCoverage.txt`, `vaccineTargetPopulation.txt` |
| Imported vaccination rates | `sigma1.txt`, `sigma2.txt`, `sigmaB.txt` |
| Fitting input | `forceFit.txt` |

The setup also sets an output path for `exitFlags.txt`, which is produced during calibration rather than read as a required input.

The Italy setup uses the **imported vaccination rates** from `sigma1.txt`, `sigma2.txt`, and `sigmaB.txt` rather than calculating these rates from coverage data during each simulation.

## Calibration and outputs

The example calls `sim.fit` to calibrate and overwrites the initial `Italy/data/kvalue.csv`.

## Visualization

The example runs `plotBothData(sim)` and `sim.fitPlot`. Additional methods, such as `sim.immunity` and `sim.compartments`, may be enabled if the corresponding plotting functions are present in `functions/`.


## Notes on reproducibility

- The Italy setup specifies a time-varying case-detection fraction and hospitalization fitting fraction, as well as Italy-specific contact patterns, vaccination rates, and length-of-stay parameters.
- Genetic algorithm estimates can differ between runs. Record the random seed, MATLAB/toolbox versions, input data, and repository commit used.

## Citation

Please cite the relevant ASIIM methodology, Italy modelling study, and the repository release or commit used.

For questions, open an issue in the [ASIIM repository](https://github.com/MarwaTuffaha/ASIIM/issues).
