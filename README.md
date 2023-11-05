# Thesis log
## Heme-Hemopexin
In the Heme-Hemopexin directory you will find three R scripts:
* Analysis - heme 1000
* Baseline method
* spatialDE pipeline - SVG

#### Significant genes
The _Significant genes_ directory contains the script with the code I used to extract significant genes (`Extract sig genes.R`) and four CSV file with the resulting genes before and after batch correction:
1. Anticorrelated_before
2. Anticorrelated_after
3. Correlated_before
4. Correlated_after

#### Supplementary 
This directory contains the following CSV files with the genes obtained by three different methods. These results were compared with each other to validate our method:
1. Down regulated - Baseline -> downreg genes using the baseline method
2. Down regulated - SVG -> downreg genes following the `spatialDE` pipeline
3. Down regulated - Spearman -> downreg genes using our method (RGDIST)
4. Up regulated - Baseline -> baseline upreg genes 
5. Up regulated - SVG -> upreg genes from `spatialDE`
6. Up regulated - Spearman -> upreg genes using RGDIST

In addition, there are the following two files containing the up- and down-regulated genes according to the authors of the dataset. These genes were used as ground truth to compare the results obtained with the three methods. 
1. heme-response down regulated genes
2. heme-response up regulated genes

## SpatialLIBD
In the [spatialLIBD](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/tree/main/spatialLIBD) directory you can find my first attempts with spatial transcriptomics data made by following the pipelines of spatialLIBD (R package v1.10.1) and SpatialExperiment (R package v1.8.1) both available on Bioconductor: _[spatialLIBD](https://bioconductor.org/packages/release/data/experiment/html/spatialLIBD.html), [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html)_.

## Wiki
In the [Wiki](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/wiki) section you can find all the updates about how to perform the analyses as well as links to interesting papers and datasets on spatial transcriptomics that I have collected throughout these months.

## Dataset
All data are available on Gene Expression Omnibus under the accession number **GSE182127**.
The dataset can also be downloaded from Zenodo at this [link](https://doi.org/10.5281/zenodo.5638720).
