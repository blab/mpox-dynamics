# monkeypox-build

## forked from [nextstrain monkeypox build](https://github.com/nextstrain/monkeypox) for use in the [monkeypox dynamics project](https://github.com/blab/monkeypox-dynamics). 

please see the nextstrain repo for documentation and usage.  

## subsampling methods for [config.yaml](config/config_hmpxv1.yaml) and [core.smk](workflow/snakemake_rules/core.smk)

Goal: retain only outbreak sequences, enrich for geographic diversity, avoid singleton countries, sample evenly through time. 

   - `min_date: 2022`
   - `--exclude-where outbreak!=hMPXV-1`
   - `--exclude-ambiguous-dates-by month` 
   - `--subsample-seed 1234567`  
   - `--sequences-per-group 5` 
   - exclude all non-B lineages
   - exclude all countries with < 3 sequences and accessions that fail downstream alignment steps 

 Note: for rooted tree, `--include-where accession=ON676708` (lineage A.1.1) in core.smk rule `exclude_bad` and rule `filter`
