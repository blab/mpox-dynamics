# Phylodynamic Analysis 

## About the Analysis

To understand mpox transmission within and between regions studied, we employ an approximate structured coalescent approach called MASCOT-GLM which also uses generalized linear models to inform estimates of effective population size and migration rates. More information can found in [Müller et al](https://academic.oup.com/ve/article/5/2/vez030/5549805?login=false)

## Running the analysis

### main analysis using GLM

1. import  `500_unmasked_equal_temporal_mascot.fasta` from [`alignments`](alignments/) into BEAUTI
2. infer tip dates 
3. non-default model parameters: 
   - SITES/site heterogeneity model: gamma with 4 sites
   - TREES/tree prior: MASCOT-GLM
   - empirical predictors can be found in [`data`](data/) which includes `region_cases_weekly_prevalence.csv`, `month_predictor_1-9.csv`, and `regions_air_travel_matrix.csv`
   - PRIORS/clock.rate: uniform prior with lower 0, upper 1.0, initial `6E-5` 
   - MCMC/length of chain: `50000000`
   - MCMC/log parameters every: `25000`
5. export xml 

All xmls used in this project can be found in [`xmls`](xmls/) and include manual xml edits to add in CoupledMCMC for faster convergence

All results for main and supplementary analyses can be found in [`results`](results/). 

Maximum clade credibility trees and effective population sizes can be analyzed [here](https://github.com/blab/mpox-dynamics/blob/main/scripts/glm_plotting_and_ne.ipynb) while percentage of cases due to introductions and Rt over time can be analyzed [here](https://github.com/blab/mpox-dynamics/blob/main/scripts/percent_cases_from_intros.ipynb)
  
XMLs and results for MASCOT-Skyline which is an approximation of the structured coalescent with a Skyline prior but without the influence of empirical predictors can be found in both `xmls/` and `results/` under the subfolder `skyline/`

