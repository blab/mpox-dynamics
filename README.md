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
   - exclude all non-B lineages
   - /add specific subsampling methods here/ 
   
4. run nextstrain workflow to produce `alignment.fasta` and `metadata.tsv` 
5. reformat `alignment.fasta` strain names -> name_location_date 

### xml preparation using BEAUTI v. 1.10.4 

1. import descriptive alignment.fasta 
2. infer tip dates 
3. non-default model parameters: 
   - SITES/site heterogeneity model: gamma
   - TREES/tree prior: coalescent bayesian skyline 
   - PRIORS/clock.rate: uniform prior with lower 0, upper 1.0, initial 6x10^-5 
   - MCMC/length of chain: 50000000
   - MCMC/log parameters every: 25000
5. export `hmpxv_skyliine.xml` 
  

### run coalescent bayesian skyline model using BEAST v. 1.10.4 




