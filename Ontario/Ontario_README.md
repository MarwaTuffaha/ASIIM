
# Ontario COVID-19 Application

This folder contains the MATLAB code and input data for applying the **Age-Structured Immunity and Infection Model (ASIIM)** to COVID-19 in Ontario, Canada. It uses the shared model and calibration functions in the repository's [`functions/`](../functions/) directory.

The example fits the model to **both reported COVID-19 cases and hospital admissions** for the pre-Omicron analysis, using a rolling-window genetic algorithm to estimate a time-varying transmission factor, $\kappa(t)$.

## Folder Structure

```text
ASIIM/
├── functions/
│   └── Shared ASIIM MATLAB functions
├── Ontario/
│   ├── example_ontario.m
│   ├── ontarioSetup.m
│   ├── data/
│   ├── Kappas from 10 calibrations/
│   └── README.md
└── README.md
```

## Requirements

- MATLAB R2021a or later (compatibility with this version has not been fully tested)
- Global Optimization Toolbox (for genetic algorithm calibration using `ga`)
- Statistics and Machine Learning Toolbox (for statistical functions such as `normrnd`)
- Parallel Computing Toolbox (if parallel optimization is enabled)

## Running the Ontario Example

1. Clone or download the [ASIIM repository](https://github.com/MarwaTuffaha/ASIIM).
2. Open MATLAB and navigate to the `Ontario` folder.
3. Run:

   ```matlab
   example_ontario
   ```

The example adds the shared functions to the MATLAB path, initializes the Ontario simulation, sets the calibration parameters, and executes:

```matlab
sim.fit
sim.simulate
sim.fitPlot
sim.immunity
```

**Each execution performs one independent calibration.** The calibration uses successive fitting windows to estimate the time-varying transmission factor.

## Input Data

The `data/` folder contains the input files required by `ontarioSetup.m`.

### Epidemiological Data
- `cases_data.csv`
- `hospital_data_unified.csv`
- `dateList.csv`

### Population and Model Parameters
- `population.txt`
- `omega.txt`
- `pax.txt`
- `kvalue.csv`

### Contact Matrices and Modifiers
- `contactModifiers.csv`
- `contactHome.csv`
- `contactSchool.csv`
- `contactWork.csv`
- `contactOtherPlaces.csv`
- `contactSF1.txt` through `contactSF3.txt`
- `contactWF1.txt` through `contactWF3.txt`
- `contactOF1.txt` through `contactOF3.txt`

### Vaccination Data
- `vaccinate.csv`
- `vaccineCoverage.csv`
- `vaccineTargetPopulation.txt`

### Calibration Inputs
- `forceFit.csv`

The file `data/kvalue.csv` provides the initial values of the time-varying transmission factor.

## Calibration

The model is calibrated to reported COVID-19 incidence and hospital admissions using a genetic algorithm applied over successive time windows.

Calibration is performed using:

```matlab
sim.fit
```

The example saves the calibrated transmission factor to:

```text
Ontario/results/kvalue_calibrated.csv
```

The original input file `data/kvalue.csv` is preserved.

## Previously Calibrated Results

The `Kappas from 10 calibrations/` folder contains ten previously obtained calibration results:

```text
kvalue1.csv
kvalue2.csv
kvalue3.csv
kvalue4.csv
kvalue5.csv
kvalue6.csv
kvalue7.csv
kvalue8.csv
kvalue9.csv
kvalue10.csv
```

These files contain previously estimated transmission factors from ten independent calibrations.

They are provided as reference results and are not automatically recalculated or loaded by the example script.

## Simulation and Visualization

After calibration, the model can be simulated using:

```matlab
sim.simulate
```

The fitted results can be visualized using:

```matlab
sim.fitPlot
```

Additional visualizations are available through the ASIIM class, provided their supporting plotting functions are included:

```matlab
sim.immunity
sim.seroprevalence
sim.compartments
```

## Notes

- Ontario-specific model configurations are defined in `ontarioSetup.m`.
- Calibration settings and execution are specified in `example_ontario.m`.
- The shared ASIIM functions are located in the repository's `functions/` directory.
- Genetic algorithm results may vary between independent calibration runs.
- For reproducibility, record the MATLAB version, toolbox versions, random seed, and repository commit.

## Citation

If you use this code, please cite the relevant ASIIM methodology and the associated Ontario COVID-19 modelling study.

Repository: [MarwaTuffaha/ASIIM](https://github.com/MarwaTuffaha/ASIIM)

## Contact

For questions or issues, please use the [GitHub Issues page](https://github.com/MarwaTuffaha/ASIIM/issues).
