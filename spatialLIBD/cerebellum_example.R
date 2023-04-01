# Using spatialLIBD with 10x Genomics public datasets. This is made by following the instructions given at the URL below
# https://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html
# For this example, I used the "Human Cerebellum: Whole Transcriptome Analysis" Visium dataset available at the link below
# https://www.10xgenomics.com/resources/datasets/human-cerebellum-whole-transcriptome-analysis-1-standard-1-2-0

library("BiocFileCache")
library("SpatialExperiment")
library("rtracklayer")
library("lobstr")
library("spatialLIBD")

## Download and save a local cache of the data provided by 10x Genomics
bfc <- BiocFileCache::BiocFileCache()
h_cereb.url <-
  paste0(
    "https://cf.10xgenomics.com/samples/spatial-exp/",
    "1.2.0/Parent_Visium_Human_Cerebellum/",
    c(
      "Parent_Visium_Human_Cerebellum_filtered_feature_bc_matrix.tar.gz",
      "Parent_Visium_Human_Cerebellum_spatial.tar.gz",
      "Parent_Visium_Human_Cerebellum_analysis.tar.gz"
    )
  )
h_cereb.data <- sapply(h_cereb.url, BiocFileCache::bfcrpath, x = bfc)

## Extract the files to a temporary location 
cereb_names <- sapply(h_cereb.data, utils::untar, exdir = file.path(tempdir(), "outs"))
## The names are the URLs, which are long and thus too wide to be shown here,
## so we shorten them to only show the file name prior to displaying the utils::untar() output status
names(cereb_names) <- basename(names(cereb_names))

## List the files we downloaded and extracted. These files are typically SpaceRanger outputs
h_cereb.dirs <- file.path(
  tempdir(), "outs",
  c("filtered_feature_bc_matrix", "spatial", "raw_feature_bc_matrix", "analysis")
)

## Import the data as a SpatialExperiment object
spe <- SpatialExperiment::read10xVisium(
  samples = tempdir(),
  sample_id = "cerebellum",
  type = "sparse", 
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

## Add some information used by spatialLIBD
spe <- add_key(spe)
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

## Initially we don't have much information about the genes. So we have to 
## read in this information from a gene annotation file: typically a gtf file. 

## For this example we download the Gencode v32 GTF file and then cache it. Depending on the version of spaceranger used we can also use 10x Genomics GTF files available at https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
gtf_cache <- BiocFileCache::bfcrpath(
  bfc,
  paste0(
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    "release_32/gencode.v32.annotation.gtf.gz"
  )
)

## Import into R (takes ~1 min)
gtf <- rtracklayer::import(gtf_cache)

## Subset to genes only
gtf <- gtf[gtf$type == "gene"]

## Remove the .x part of the gene IDs
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)

## Set the names to be the gene IDs
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(spe), gtf$gene_id)
table(is.na(match_genes))

## Drop the few genes for which we don't have information
spe <- spe[!is.na(match_genes), ]
match_genes <- match_genes[!is.na(match_genes)]

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]


## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

## Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
## Number of genes with no counts
#  length(no_expr)

## Compute the percent of genes with no counts
#> length(no_expr) / nrow(spe) * 100
spe <- spe[-no_expr, , drop = FALSE]

## Remove spots without counts
summary(spe$sum_umi)
## If we had spots with no counts, we would remove them
if (any(spe$sum_umi == 0)) {
  spots_no_counts <- which(spe$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

## Add a variable for saving the manual annotations
spe$ManualAnnotation <- "NA"

## The logNormCounts() function from scuttle will compute a log-transformed normalized 
## expression matrix and store it as another assay ("logcounts"). 
spe <- scuttle::logNormCounts(spe)

# txt file containing the genes ID stored in rowDAta(spe)
write.table(rowData(spe)$gene_search, "C:/Users/mason/OneDrive/Documenti/genes_ID.txt" )


# "Many genes, including En1, En2, Pax2, Wnt7b, and some of the ephrins and their receptors, 
# show characteristic patterns of spatial expression in the cerebellum, but only En2 has been 
# studied specifically for its role in compartmentalization." https://www.nature.com/articles/35081558#:~:text=Many%20genes%2C%20including%20En1%2C%20En2,for%20its%20role%20in%20compartmentalization

## Example visualizations of En2
spatialLIBD::vis_gene(
  spe = spe,
  sampleid = "cerebellum",
  geneid = "EN2; ENSG00000164778",
  assayname = "logcounts",
  viridis = FALSE
)


