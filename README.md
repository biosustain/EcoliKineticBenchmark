# E. coli kinetic model benchmark

## Overview
Repository contains:
1) Mathematical models in SBML or Matlab format that can simulate E. coli metabolism - `data/models/*`.
2) Results of such simulations - `data/simulation_results`.
3) Experimental data from publications with similar experimental design - `data/datasets/*`.
4) Set of notebooks to orchestrate simulations and analysis of results - `notebooks/`.
5) Additional files that are dictionaries to convert IDs from model specific ones to the BiGG IDs.

## Instructions
Most of the tasks should be started from the notebooks titled `Run Sensitivity experiments` and `Run KO experiment`.
COBRA simulation is currently being done from the notebook `COBRA Simulation`.

## Requirements
`python > 3.5`
`tellurium > 2.0`
`pandas > 0.20`
`numpy > 1.14`
`altair > 2.1`
`cobrapy > 0.14`