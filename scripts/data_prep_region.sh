# multiple sequence alignment file prep for beauti
# after completing nextstrain workflow through 'mask' steps,
# make sure there are no duplicate strain names prior to running
# usage:
# ./data_prep.sh


mkdir ../monkeypox-build/results/beauti
cp ../monkeypox-build/data/metadata_region_1000.tsv ../monkeypox-build/results/beauti/
cp ../monkeypox-build/results/hmpxv1/masked.fasta ../monkeypox-build/results/beauti/
cd ../monkeypox-build/results/beauti

#make new strain name col with format name_region_date
awk -F"\t" 'OFS="\t" {$1=$6"_"$4"_"$3; print}' metadata_region_1000.tsv > meta.tsv

#remove spaces from strain names
awk -F"\t" 'OFS="\t" {gsub(/[[:blank:]]/, "",$1); print}' meta.tsv > tmp && mv tmp meta.tsv

#rename column as 'strain'
awk -F"\t" 'OFS="\t" {sub(/strain_region_date/,"strain",$1); print}' meta.tsv > tmp && mv tmp meta.tsv

#make key,value file kv.txt with 1st col 'accessions', 2nd col 'strain'
awk 'NR==1{OFS="\t";save=$2;print $2,$1}NR>1{print $2,$1,save}' meta.tsv > kv.txt
awk '!($3="")' kv.txt

#replace old strain names in alignment .fasta  with new strain names
cat masked.fasta | seqkit replace --ignore-case --kv-file "kv.txt" --pattern "^(\w+)" --replacement "{kv}" > align.fasta

#remove the outgroup for beast analyses
cat align.fasta | awk '/>MPXV_USA_2021_MD_NorthAmerica_2021-11-25/ {getline; while(!/>/) {getline}} 1' > fixed_region_1000.fasta

