## --- Making a single spe object containing all samples ---
# To work with multiple samples we need to create a spe object for each one 
# and then combine them into a single spe object using the cbind() function

# Control: mouse brain without any treatment
spe_c <- SpatialExperiment::read10xVisium(
  samples = "control",
  sample_id = "control",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# heme_0030: 0.30 nmol heme-albumin in 10 µL
spe_30 <- SpatialExperiment::read10xVisium(
  samples = "heme_0030",
  sample_id = "heme_0030",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# heme_0125: 1.25 nmol heme-albumin in 10 µL
spe_125 <- SpatialExperiment::read10xVisium(
  samples = "heme_0125",
  sample_id = "heme_0125",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# heme_0500: 5.0 nmol heme-albumin in 10 µL
spe_500 <- SpatialExperiment::read10xVisium(
  samples = "heme_0500",
  sample_id = "heme_0500",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# heme_1000: 10 nmol heme-albumin in 10 µL
spe_1000 <- SpatialExperiment::read10xVisium(
  samples = "heme_1000",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# hemeHpx_1000: 10 nmol heme-hemopexin in 10 µL
spe_hpx <- SpatialExperiment::read10xVisium(
  samples = "hemeHpx_1000",
  sample_id = "hemeHpx_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# sham: 10 µL saline
spe_s <- SpatialExperiment::read10xVisium(
  samples = "sham",
  sample_id = "sham",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# combine object
spe <- cbind(spe_c, spe_30, spe_125, spe_500, spe_1000, spe_hpx, spe_s)

# modify spe to spatialLIBD
spe <- add_key(spe)
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
spe <- scuttle::logNormCounts(spe)

gtf <-                                           
  rtracklayer::import(
    "gencode.vM25.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    
names(gtf) <- gtf$gene_id                        

match_genes <- match(rownames(spe), gtf$gene_id) 
table(is.na(match_genes))

spe <- spe[!is.na(match_genes), ]                
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  
rowRanges(spe) <- gtf[match_genes]               

rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

no_expr <- which(rowSums(counts(spe)) == 0)           
length(no_expr) / nrow(spe) * 100                     
spe <- spe[-no_expr, , drop = FALSE]

summary(spe$sum_umi)

if (any(spe$sum_umi == 0)) {
  spots_no_counts <- which(spe$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

spe$ManualAnnotation <- "NA"                         
check_spe(spe)

# Hmox1; ENSMUSG00000005413
p_list <-
  vis_grid_gene(
    spe[, spe$sample_id %in% c("control", "sham", "heme_1000")],
    geneid = "Hmox1; ENSMUSG00000005413",
    assayname = "counts",
    return_plots = TRUE,
    spatial = FALSE,
  )

cowplot::plot_grid(plotlist = p_list, ncol = 3)

