# E. coli kinetic model benchmark
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3610129.svg)](https://doi.org/10.5281/zenodo.3610129)
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

## References

Please cite original papers if you use any of their data or repackaged versions from this repository in your projects.

The data that was used:
1) Toya et al., “13C-Metabolic Flux Analysis for Batch Culture of Escherichia Coli and Its Pyk and Pgi Gene Knockout Mutants Based on Mass Isotopomer Distribution of Intracellular Metabolites.”
2) Nicolas et al., “Response of the Central Metabolism of Escherichia Coli to Modified Expression of the Gene Encoding the Glucose-6-Phosphate Dehydrogenase.”
3) Yao et al., “Catabolic Regulation Analysis of Escherichia Coli and Its Crp, Mlc, MgsA, Pgi and PtsG Mutants.”
4) Usui et al., “Investigating the Effects of Perturbations to Pgi and Eno Gene Expression on Central Carbon Metabolism in Escherichia Coli Using (13)C Metabolic Flux Analysis.”
5) Long and Antoniewicz, “Metabolic Flux Responses to Deletion of 20 Core Enzymes Reveal Flexibility and Limits of E. Coli Metabolism.”

The models:
1) Chassagnole et al., “Dynamic Modeling of the Central Carbon Metabolism of Escherichia Coli.”
2) Khodayari and Maranas, “A Genome-Scale Escherichia Coli Kinetic Metabolic Model k-Ecoli457 Satisfying Flux Data for Multiple Mutant Strains.”
3) Millard, Smallbone, and Mendes, “Metabolic Regulation Is Sufficient for Global and Robust Coordination of Glucose Uptake, Catabolism, Energy Production and Growth in Escherichia Coli.”
4) Kurata and Sugimoto, “Improved Kinetic Model of Escherichia Coli Central Carbon Metabolism in Batch and Continuous Cultures.”