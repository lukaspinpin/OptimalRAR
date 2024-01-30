# OptimalRAR: Simulation of Response-Adaptive Randomization Designs

This repository contains the source code used for the simulations presented in the book chapter "Response-Adaptive Randomization Designs based on Optimal Allocation Proportions". The codebase includes various scripts to simulate different response-adaptive randomization scenarios and analyze the results.

## Contents

- `Simulation.R`: This is the main script of the project. By default, it is set up to replicate the results for Table 3 from the book chapter. However, it is highly configurable, allowing users to simulate different scenarios by adjusting parameters such as sample size, burn-in period, and endpoints.

- `Data_Generator.R`: This script generates data based on user-specified parameters. It forms the backbone of the simulation process, providing the necessary data for analysis.

- `WaldTest.R`: Implements the Wald Test for the simple mean difference across various endpoints. It includes specific functions for binary endpoints and the six measures discussed in the book chapter. The test statistic is also corrected for extreme cases.

- `ER.R`: Script for Equal Randomization (ER). To implement restricted randomization and ensure equal sample sizes per arm, set the burn-in period to half of the sample size (`n/2`).

- `SQE.R`: Implements SMLE targeting either Neyman or minF proportions.

- `DBCD.R`: DBCD also targeting Neyman or minF proportions.

- `ERARDE.R`: ERADE targeting Neyman or minF proportions.

## Usage

The scripts are designed to be modular and easily customizable. Users can modify the parameters in `Simulation.R` to explore different adaptive randomization strategies and outcomes. For detailed instructions on how to set up and run the simulations, please refer to the comments within each script.

## Erros 
If you believe you have identified any errors, please reach out to Lukas Pin at lukas.pin@mrc-bsu.cam.ac.uk.
