Anaconda with python 3.6 recommended, with environment sourced
R 3.5.1

Dataset: TCGA 450K

# Download pymethylprocess
Run the following commands:

```
pip install git+https://github.com/bodono/scs-python.git@bb45c69ce57b1fbb5ab23e02b30549a7e0b801e3 git+https://github.com/jlevy44/hypopt.git@af59fbed732f5377cda73fdf42f3d4981c2be3ce
pip install pymethylprocess && pymethyl-install_r_dependencies
```

# Running the preprocessing pipeline:


Download TCGA files:
```
nohup pymethyl-preprocess download_tcga
pymethyl-preprocess download_clinical
```
Create and format sample sheets:
```
pymethyl-preprocess create_sample_sheet -is ./tcga_idats/clinical_info.csv -s tcga -i tcga_idats/ -os tcga_idats/samplesheet.csv
mkdir backup_clinical && mv ./tcga_idats/clinical_info.csv backup_clinical
```
Filter out subtypes with counts below a threshold:
```
pymethyl-preprocess get_categorical_distribution -k disease
pymethyl-preprocess remove_diseases -is ./tcga_idats/samplesheet.csv -os ./tcga_idats/samplesheet_exclusion.csv -l 11 -d
```
Format phenosheets:
```
pymethyl-preprocess meffil_encode -is ./tcga_idats/samplesheet_exclusion.csv -os ./tcga_idats/final_samplesheet.csv
```
Split phenosheet by disease subtype ad run preprocessing for each subtype in parallel (combine arrays at end):
```
pymethyl-preprocess split_preprocess_input_by_subtype -i ../pancancer/final_samplesheet.csv -d
nohup time pymethyl-preprocess batch_deploy_preprocess -n 6 -c 5 -r -m &
pymethyl-preprocess combine_methylation_arrays -d ./preprocess_outputs/ -e OV
```
Remove non-autosomal CpGs, SNPs and report missingness:
```
pymethyl-utils print_number_sex_cpgs -i combined_outputs/methyl_array.pkl
pymethyl-utils remove_sex -i combined_outputs/methyl_array.pkl
pymethyl-utils remove_snps -i autosomal/methyl_array.pkl
pymethyl-preprocess na_report -i autosomal/methyl_array.pkl -o na_report/
```
Run imputation pipeline with feature selection using mean absolute deviation:
```
nohup pymethyl-preprocess imputation_pipeline -i ./no_snp/methyl_array.pkl -o final_preprocessed/methyl_array.pkl -n 200000 -ss -s sklearn -m Mean -k 20 -d -st 0.05 -ct 0.05 &
```
Generate visualizations:
```
mkdir visualizations
nohup pymethyl-visualize transform_plot -o visualizations/umap_embed.html -c disease_only -nn 8 &
```
Split dataset into training, testing, and validation sets:
```
pymethyl-utils train_test_val_split -tp .8 -vp .125 -cat
pymethyl-utils counts -i train_val_test_sets
```
