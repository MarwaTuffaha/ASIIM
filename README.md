
# ASIIM-COVID19

MATLAB source code for the **Age-Structured Immunity and Infection Model (ASIIM)**, with applications to COVID-19 modelling in different jurisdictions.

## Overview

This repository contains the source code used to implement ASIIM and apply it to COVID-19 epidemiological modelling across different jurisdictions.

ASIIM uses a common set of MATLAB functions that support different model configurations. Each jurisdiction has its own scripts and inputs that use these shared functions.

## Repository Structure

```text
ASIIM-COVID19/
├── functions/       # Shared MATLAB functions
├── jurisdiction1/   # scripts and configurations specific to jurisdiction 1
├── jurisdiction2/   # scripts and configurations specific to jurisdiction 2
├── ....
└── README.md
```

### Shared Functions

The `functions/` directory contains the common MATLAB functions implementing ASIIM. These functions are used across all jurisdiction-specific applications.

### Jurisdiction-Specific Code

Each jurisdiction directory contains scripts and configurations specific to its application:

All applications use the shared ASIIM functions rather than maintaining separate implementations of the model.

## Requirements

- MATLAB
- Any additional MATLAB toolboxes required by the relevant scripts.

## Getting Started

1. Clone or download this repository.
2. Open MATLAB and navigate to the repository directory.
3. Add the shared functions to the MATLAB path:

   ```matlab
   addpath(genpath('functions'));
   ```

4. Navigate to the relevant jurisdiction directory.
5. Run the appropriate scripts using the provided configurations and inputs.

Execution steps and required inputs may vary between jurisdictions.

## Reproducibility

The repository provides a common implementation of ASIIM, with jurisdiction-specific configurations and analysis scripts.

For reproducibility, users are encouraged to reference the relevant publication and the repository version or commit used.

## Related Publications

Publications associated with the Ontario, Italy, and United States applications will be listed here as they become available.

## Citation

If you use this source code, please cite the relevant ASIIM publication and the corresponding jurisdiction-specific study.

## Contact

For questions or issues related to the source code, please open an issue in this GitHub repository.
