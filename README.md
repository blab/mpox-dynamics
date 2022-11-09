# monkeypox-dynamics:

## phylodynamic inference of hmpxv1 effective population size during 2022 outbreak

## analysis

### data preparation 

1. clone [nexstrain monkeypox build](https://github.com/nextstrain/monkeypox) |  [my fork](https://github.com/nmmahmed/monkeypox-build) 
2. fetch data from nextstrain.org 
3. modify filter rule in config file 
   - `min_date: 2022`
   - `--exclude-where outbreak!=hMPXV-1`
   - `--exclude-ambiguous-dates-by month` 
   - `--subsample-seed 1234567`  
   - exclude all non-B lineages
   - /add my specific subsampling methods here/
   
4. run nextstrain workflow to produce `alignment.fasta` and `metadata.tsv` 
5. [mask](/scripts/masking.ipynb) invariant sites from alignment file 
   - modify mask rule in config file to accept new BED file
6. re-run nextstrain workflow using identical `subsample-seed` 
5. [reformat](/scripts/data_prep.sh) `alignment.fasta` strain names as name_location_date

### xml preparation using BEAUTI v. 1.10.4 

1. import descriptive `alignment.fasta` 
2. infer tip dates 
3. non-default model parameters: 
   - SITES/site heterogeneity model: gamma
   - TREES/tree prior: coalescent bayesian skyline 
   - PRIORS/clock.rate: uniform prior with lower 0, upper 1.0, initial `6E-5` 
   - MCMC/length of chain: `50000000`
   - MCMC/log parameters every: `25000`
5. export `hmpxv_skyliine.xml` 
6. delete 'trait' block containing 'strictClockBranchRates'
  
### run coalescent bayesian skyline model using BEAST (v. 1.10.4) 

### estimate TMRCA and plot skyline using TRACER (v. 1.7.2) 

### run coalescent logistic growth model using BEAST (v. 1.10.4) 

XML preparation 
   - SITES/site heterogeneity model: gamma
   - TREES/tree prior: coalescent logistic growth model 
   - PRIORS/clock.rate: uniform prior with lower 0, upper 1.0, initial `6E-5` 
   - MCMC/length of chain: `50000000`
   - MCMC/log parameters every: `25000`
Export `hmpxv_skyline_logistic.xml`. Delete 'trait' block containing 'strictClockBranchRates'; reformat 'gammaPrior' block as 'uniformPrior'. 

### 








