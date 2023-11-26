# Discete Trait Analysis 

## About the Analysis

To understand global spread of mpox, we employ phylogeographic models to estimate the rates of migration between our five studied regions.

## Running the analysis

### main analysis using Skygrid prior

1. import  `main_masked_prevalence_subsampling_alignment.fasta` from [`alignments`](alignments/) into BEAUTI (v. 1.10.4)
2. infer tip dates 
3. non-default model parameters: 
   - SITES/site heterogeneity model: gamma with 4 sites
   - TREES/tree prior: coalescent bayesian Skygrid with 1.5 years before latest tip and 75 change points to allow for the effective population size to change every two weeks 
   - PRIORS/clock.rate: uniform prior with lower 0, upper 1.0, initial `6E-5` 
   - MCMC/length of chain: `50000000`
   - MCMC/log parameters every: `25000`
5. export xml 

All xmls used in this project can be found in [`xmls`](xmls/). 

All results for main and supplementary analyses can be found in [`results`](results/)

Maximum clade credibility trees can be analyzed [here](https://github.com/blab/mpox-dynamics/blob/main/scripts/plotting_dta_tree.ipynb) while imports, exports, and persistence times can be analyzed [here](https://github.com/blab/mpox-dynamics/blob/main/scripts/estimating_introductions_exportations.ipynb)
  
