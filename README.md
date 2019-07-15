# E. coli kinetic model benchmark

## Overview
Repository contains:
1) Mathematical models in SBML or Matlab format that can simulate E. coli metabolism - `data/models/*`.
2) Results of such simulations - `data/simulation_results`.
3) Experimental data from publications with similar experimental design - `data/datasets/*`.
4) Set of notebooks to orchestrate simulations and analysis of results - `notebooks/`.
5) Additional files that are dictionaries to convert IDs from model specific ones to the BiGG IDs.

## Additional links
Models overview and scope - [Here](https://docs.google.com/spreadsheets/d/1AabFcblykcDXuCHVakR9OE9DLRcwop4-ajbPAx--nSI/edit?usp=sharing)
Data collection - [Here](https://docs.google.com/spreadsheets/d/1P_F2AuiC-Hcu-Z4hUZIyO9lj60QOr3mtdfIW9VuNxqs/edit?usp=sharing)

## Instructions
Most of the tasks should be started from the notebooks titled `Run Sensitivity experiments` and `Run KO experiment`.
COBRA simulation is currently being done from the notebook `COBRA Simulation`.

## Requirements
```
python > 3.6
tellurium > 2.0
pandas > 0.24
numpy > 1.16
altair > 3.0
cobrapy > 0.15
```