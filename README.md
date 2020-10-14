# precise-db

PRECISE is an RNA-seq compendium for _Escherichia coli_ from the Systems Biology Research group at UC San Diego.

[![DOI](https://zenodo.org/badge/187104662.svg)](https://zenodo.org/badge/latestdoi/187104662)

## Data
The following data files are available in the `data` folder:
* log_tpm_full.csv: Expression levels for all genes in E. coli
* log_tpm.csv: Expression levels for 3,923 genes in E. coli (noisy genes have been removed)
* log_tpm_norm.csv: `log_tpm.csv` centered to reference condition (WT on glucose M9 media)
* metadata.csv: Experimental metadata (e.g. strain descriptions, carbon source etc.) for all 278 conditions in PRECISE
* gene_info.csv: Descriptive characteristics of genes, including location, operon, and COG group
* TRN.csv: Known regulator-gene interactions from RegulonDB 10.0
* S.csv: Gene coefficients for each iModulon
* A.csv: Condition-specific activities for each iModulon
* curated_enrichments.csv: Detailed information on iModulons and their linked regulator(s)
* imodulon_gene_names.txt: List of gene names in each iModulon
* imodulon_gene_bnumbers.txt: List of genes (as b-numbers) in each iModulon

## Scripts
A conda environment for this code has been provided [here](https://github.com/SBRG/precise-db/blob/master/environment.yml)  
To generate robust independent components for a dataset, execute the `run_ica.sh` script:  
`run_ica <filename.csv>`  
where `<filename.csv>` is a comma-separated file of gene expression. Data must be centered using a reference condition (See `data/log_tpm_norm.csv` for an example)
Additional options are included as flags. Decreasing tolerance (e.g. `-t 1e-3`) will reduce runtime, but will also reduce the final number of independent components.

## Notebooks
The Jupyter notebook `exploratory_analysis.ipynb` walks users through the data files and includes a few small functions for interrogating iModulons.

## Requirements:
Python 3.6 or greater  
Conda environment specifications are listed in [environment.yml](https://github.com/SBRG/precise-db/blob/master/environment.yml)   
Versions of `scikit-learn` above `0.20.3` cause an error when performing ICA.
