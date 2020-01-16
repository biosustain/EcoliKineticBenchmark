# E. coli kinetic model benchmark

## Overview
Repository contains:
1) Mathematical models in SBML or Matlab format that can simulate E. coli metabolism - `data/models/*`.
2) Results of such simulations - `data/simulation_results`.
3) Experimental data from publications with similar experimental design - `data/datasets/*`.
4) Set of notebooks to orchestrate simulations and analysis of results - `notebooks/`.
5) Additional files that are dictionaries to convert IDs from model specific ones to the BiGG IDs `data/`.


## Instructions
Most of the tasks should be started from the notebooks titled `Run Sensitivity experiments` and `Run KO experiment`.

COBRA simulation is currently being done from the notebook `COBRA Simulation`.

After simulations are complete use notebooks `notebooks/Analyze *` to generate required visualizations.


## Requirements
```
python > 3.6
tellurium > 2.1
pandas > 0.24
numpy > 1.16
xarray > 0.12.0
altair > 3.2.0
cobrapy > 0.15
escher = 1.6.0
```