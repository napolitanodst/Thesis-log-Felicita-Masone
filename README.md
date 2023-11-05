# Thesis log

## Heme-Hemopexin
In the 'Heme-Hemopexin' directory, you will find three R scripts and two subfolders:
* Analysis - heme 1000.R
* Baseline method.R
* spatialDE pipeline - SVG.R
* Significant genes
* Supplementary

#### Significant genes
The 'Significant Genes' directory contains the script (`Extract sig genes.R`) that I used to extract significant genes, as well as four CSV files with the resulting genes before and after batch correction.:
1. Anticorrelated_before
2. Anticorrelated_after
3. Correlated_before
4. Correlated_after

#### Supplementary 
This directory contains the following CSV files with genes obtained through three different methods. These results were compared to validate our method:
1. Downregulated - Baseline: Genes downregulated using the baseline method.
2. Downregulated - SVG: Genes downregulated following the `spatialDE` pipeline.
3. Downregulated - Spearman: Genes downregulated using our method (RGDIST).
4. Upregulated - Baseline: Genes upregulated using the baseline method.
5. Upregulated - SVG: Genes upregulated from `spatialDE`.
6. Upregulated - Spearman: Genes upregulated using RGDIST.

In addition, there are two files containing genes identified as up- and down-regulated by the authors of the dataset. These genes were used as ground truth to compare the results obtained with the three methods:
1. Heme-Response Downregulated Genes
2. Heme-Response Upregulated Genes

## SpatialLIBD
In the [spatialLIBD](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/tree/main/spatialLIBD) directory you can find my first attempts with spatial transcriptomics data made by following the pipelines of spatialLIBD (R package v1.10.1) and SpatialExperiment (R package v1.8.1) both available on Bioconductor: _[spatialLIBD](https://bioconductor.org/packages/release/data/experiment/html/spatialLIBD.html), [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html)_.

## Wiki
In the [Wiki](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/wiki) section you can find all the updates about how to perform the analyses as well as links to interesting papers and datasets on spatial transcriptomics that I have collected throughout these months.

## Dataset
All data are available on Gene Expression Omnibus under the accession number **GSE182127**.
The dataset can also be downloaded from Zenodo at this [link](https://doi.org/10.5281/zenodo.5638720).

## [Identification of treatment-responsive genes in spatial transcriptomics data by leveraging injection site information](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/blob/main/Identification%20of%20treatment-responsive%20genes%20in%20spatial%20transcriptomics%20data%20by%20leveraging%20injection%20site%20information.pdf)
This is the short paper that Professor Francesco Napolitano and I presented at the CIBB conference on September 8, 2023, in Padova. The article pertains to the RGDIST method developed during my thesis work. Currently, we are actively working to enhance both the statistical and computational aspects of the method.  

## [LaTeX thesis template](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/blob/main/Template%20tesi%20LaTeX.zip)
In this zipped folder, you will find my LaTeX thesis template. To use it, simply download the folder and upload it to Overleaf.

## [Full text thesis (it)](https://github.com/napolitanodst/Thesis-log-Felicita-Masone/blob/main/Tesi%20completa.zip)
