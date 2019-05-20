# precise-db

PRECISE is an RNA-seq compendium for _Escherichia coli_ from the Systems Biology Research group at UC San Diego.

## Data
The following data files are available in the `data` folder:
* log_tpm_full.csv: Expression levels for all genes in E. coli
* log_tpm.csv: Expression levels for 3,923 genes in E. coli (noisy genes have been removed)
* metadata.csv: Experimental metadata (e.g. strain descriptions, carbon source etc.) for all 278 conditions in PRECISE
* gene_info.csv: Descriptive characteristics of genes, including location, operon, and COG group
* TRN.csv: Known regulator-gene interactions from RegulonDB 10.0
* S.csv: Gene coefficients for each i-modulon
* A.csv: Condition-specific activities for each i-modulon
* curated_enrichments.csv: Detailed information on i-modulons and their linked regulator(s)

The Jupyter notebook `exploratory_analysis.ipynb` walks users through these data files and includes a few small functions for interrogating i-modulons.

Requirements:
Python 3.6 or greater

pandas>=0.24.2  
matplotlib>=3.0.3  
scipy>=1.2.1  
