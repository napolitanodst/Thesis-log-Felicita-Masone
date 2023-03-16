# Create a spe object with the Spatial transcriptome data from coronal mouse brain sections after striatal injection of heme and heme-hemopexin. 

library("BiocFileCache")
library("SpatialExperiment")
library("rtracklayer")
library("lobstr")
library("spatialLIBD")
library("readr") 

## 1--- IMPORT THE DATA AS A SpatialExperiment OBJECT

# Control: mouse brain without any treatment
spe <- SpatialExperiment::read10xVisium(
  samples = "control",
  sample_id = "control",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

## 2--- MODIFY spe FOR spatialLIBD

# Add some information used by spatialLIBD
spe <- add_key(spe)
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
spe <- scuttle::logNormCounts(spe)


## 3--- ADD GENE ANNOTATION INFO (reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/)

gtf <-                                           # Read in the gene information from the annotation GTF file (from genecode)
  rtracklayer::import(
    "gencode.vM25.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   # Subset to genes only
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    # Remove the .x part of the gene IDs
names(gtf) <- gtf$gene_id                        # Set the names to be the gene IDs

match_genes <- match(rownames(spe), gtf$gene_id) # Match the genes
table(is.na(match_genes))

spe <- spe[!is.na(match_genes), ]                # Drop the few genes for which we don't have information
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  # Keep only some columns from the gtf

rowRanges(spe) <- gtf[match_genes]               # Add the gene info to our SPE object

# Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)
# Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

## 4--- FILTER THE spe OBJECT

no_expr <- which(rowSums(counts(spe)) == 0)           # Remove genes with no data (Number of genes with no counts: length(no_expr))
length(no_expr) / nrow(spe) * 100                     # Compute the percent of genes with no counts
spe <- spe[-no_expr, , drop = FALSE]

summary(spe$sum_umi)
# If we had spots with no counts, we would remove them
if (any(spe$sum_umi == 0)) {
  spots_no_counts <- which(spe$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

## 5--- check spe

spe$ManualAnnotation <- "NA"                                                # Add a variable for saving the manual annotations

control_anatomy <- read_csv("control/control_anatomy.csv")                  # Import anatomical information and add them to spe
spe_c$anatomy <- control_anatomy$anatomy

control_injection_site <- read_csv("control/control_injection_site.csv")    # Import the injection site data and add them to spe     
spe_c$inj_site <- control_injection_site$injection_site

check_spe(spe)

## 6--- Explore the data

# Example of spatial gene visualization of Mt2
vis_gene(
  spe = spe,
  sampleid = "control",
  geneid = "Mt2; ENSMUSG00000031762", 
  spatial = TRUE,
  assayname = "counts",
  viridis = T,
  image_id = "lowres",
  alpha = 1,
  point_size = 2,
  ... = ""
)

# Example of spatial gene visualization of Cd63
vis_gene(
  spe = spe,
  sampleid = "control",
  geneid = "Cd63; ENSMUSG00000025351", 
  spatial = TRUE,
  assayname = "counts",
  viridis = T,
  image_id = "lowres",
  alpha = 1,
  point_size = 2,
  ... = ""
)


