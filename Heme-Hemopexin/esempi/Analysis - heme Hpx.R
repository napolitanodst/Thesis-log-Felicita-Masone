# SAMPLE: hemeHpx_1000 

# Import data as SpatialExperiment object ----
library(SpatialExperiment)
library(rtracklayer)
library(lobstr)

spe_hH <- SpatialExperiment::read10xVisium(
  samples = "hemeHpx_1000",
  sample_id = "hemeHpx_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# continuous variable containing total number of counts for each sample prior to filtering any genes
spe_hH$sum_umi <- colSums(counts(spe_hH))
# continuous variable containing the number of genes that have at least 1 count
spe_hH$sum_gene <- colSums(counts(spe_hH) > 0)

# gene annotation - reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/
gtf <-                                           
  rtracklayer::import(
    "gencode.vM23.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    
names(gtf) <- gtf$gene_id                        

match_genes <- match(rownames(spe_hH), gtf$gene_id) 
table(is.na(match_genes))
spe_hH <- spe_hH[!is.na(match_genes), ]                
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  
rowRanges(spe_hH) <- gtf[match_genes]               

rowData(spe_hH)$gene_search <- paste0(
  rowData(spe_hH)$gene_name, "; ", rowData(spe_hH)$gene_id
)

# Anatomy info ----
library(readr)

hemeHpx_1000_anatomy <- read_csv("hemeHpx_1000/hemeHpx_1000_anatomy.csv")
spe_hH$anatomy <- hemeHpx_1000_anatomy$anatomy

# Distances from injection site ----

tissue_positions_list_hemeHpx_1000 <- read_csv("hemeHpx_1000/outs/spatial/tissue_positions_list.csv", col_names = FALSE)
names(tissue_positions_list_hemeHpx_1000) <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
tissue_positions_list_hemeHpx_1000 <- subset(tissue_positions_list_hemeHpx_1000, in_tissue == 1)
tissue_positions_list_hemeHpx_1000 <- tissue_positions_list_hemeHpx_1000[,-c(5,6)]
hemeHpx_1000_positions <- merge(tissue_positions_list_hemeHpx_1000, hemeHpx_1000_anatomy, by = "Barcode")
remove(tissue_positions_list_hemeHpx_1000)

# calcolo la distanza utilizzando le coordinate dello spot centrale del sito di iniezione 
# TGGCAAACTAAATTAC-1 x = 50 ; y = 78  
hemeHpx_1000_positions$distance <- NA
for(i in 1:nrow(hemeHpx_1000_positions)){
  hemeHpx_1000_positions[i,6] <- sqrt((hemeHpx_1000_positions[i,3] - 50)^2 + (hemeHpx_1000_positions[i,4] - 78)^2)
}

colData(spe_hH)$inj_site_distance <- hemeHpx_1000_positions$distance

# Quality Control and data filtering ----
library(scater)
library(ggspavis)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_hH)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe_hH <- addPerCellQC(spe_hH, subsets = list(mito = is_mito))

# histogram of numbers of expressed genes
hist(colData(spe_hH)$detected, breaks = 20)

# Remove mitochondrial genes
spe_hH <- spe_hH[!is_mito, ]

# remove not expressed genes
no_expr <- which(rowSums(counts(spe_hH)) == 0)           
length(no_expr) / nrow(spe_hH) * 100                     
spe_hH <- spe_hH[-no_expr, , drop = FALSE]

summary(spe_hH$sum)

if (any(spe_hH$sum == 0)) {
  spots_no_counts <- which(spe_hH$sum == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe_hH) * 100)
  spe_hH <- spe_hH[, -spots_no_counts, drop = FALSE]
}

# Check the number of spots in which a gene is expressed and removal of genes expressed in less than 3 spots

hH_counts <- assays(spe_hH)$counts
hH_counts <- as.matrix(hH_counts)

no_rel <- c()
n <- 0

for(i in 1:nrow(hH_counts)){
  for(j in 1:ncol(hH_counts)){
    if(hH_counts[i,j] != 0){
      n = n + 1
    }
  }
  if(n <= 2){
    no_rel = c(no_rel, i)
  }
  n = 0
}

remove(hH_counts)
spe_hH <- spe_hH[-no_rel, , drop = FALSE]

# Normalization
spe_hH <- scuttle::logNormCounts(spe_hH)

# Batch Correction on anatomy cluster ---- 
library(sva)

hH_counts <- assays(spe_hH)$counts
hH_counts <- as.matrix(hH_counts)

batch <- c()
for(i in 1:nrow(hemeHpx_1000_positions)){
  if(is.na(hemeHpx_1000_positions[i,5])){
    batch[i] <- NA
  }else{
    if(hemeHpx_1000_positions[i,5] == "caudate_putamen"){
      batch[i] <- 1
    }else{
      if(hemeHpx_1000_positions[i,5] == "cortex"){
        batch[i] <- 2
      }else{
        if(hemeHpx_1000_positions[i,5] == "thalamus"){
          batch[i] <- 3
        }else{
          if(hemeHpx_1000_positions[i,5] == "globus_pallidus"){
            batch[i] <- 4
          }else{
            if(hemeHpx_1000_positions[i,5] == "plexus"){
              batch[i] <- 5
            }else{
              if(hemeHpx_1000_positions[i,5] == "hypothalamus"){
                batch[i] <- 6
              }else{
                if(hemeHpx_1000_positions[i,5] == "corpus_callosum"){
                  batch[i] <- 7
                }
              }
            }
          }
        }
      }
    }
  }
}
hemeHpx_1000_positions$batch <- batch

# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(batch)))  
batch_noNA <- batch[!is.na(batch)]
hH_counts_noNA <- hH_counts[,-NA_spot] 

spe_hH_noNA <- spe_hH[,-NA_spot]

BCmatrix_noNA_c <- ComBat_seq(counts = hH_counts_noNA,
                              batch = batch_noNA,
                              group = NULL,
                              covar_mod = NULL,
                              full_mod = FALSE)

assays(spe_hH_noNA)$BC <- BCmatrix_noNA_hH

# Check the number of spots in which a gene is expressed and removal of genes expressed in less than 3 spots

no_relNA <- c()
n <- 0

for(i in 1:nrow(BCmatrix_noNA_hH)){
  for(j in 1:ncol(BCmatrix_noNA_hH)){
    if(BCmatrix_noNA_hH[i,j] != 0){
      n = n + 1
    }
  }
  if(n <= 2){
    no_relNA = c(no_relNA, i)
  }
  n = 0
}

spe_hH_noNA <- spe_hH_noNA[-no_relNA, , drop = FALSE]

spe_hH_noNA <- logNormCounts(x = spe_hH_noNA,
                                assay.type = "BC",
                                name = "BClogcounts")

# Reduced dimension: UMAP ----
library(scater)

# No corrected data UMAP
spe_hH <- runUMAP(spe_hH, exprs_values = "logcounts")
reducedDimNames(spe_hH)
head(reducedDim(spe_hH))
plotUMAP(spe_hH, colour_by = "anatomy")

# Batch Corrected data UMAP
spe_hH_noNA <- runUMAP(spe_hH_noNA, exprs_values = "BC")
reducedDimNames(spe_hH_noNA)
head(reducedDim(spe_hH_noNA))
plotUMAP(spe_hH_noNA, colour_by = "anatomy")

# Correlation and p.values ----

# BEFORE batch effect removal
hH_logcounts <- assays(spe_hH)$logcounts
hH_logcounts <- as.matrix(hH_logcounts)
distance_hH <- spe_hH$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_hH <- c()
p_value_hH <- c()

for(i in 1:nrow(hH_logcounts)) {
  zero <- which(hH_logcounts[i,] == 0)
  cor_hH[i] <- cor(distance_hH[-zero], hH_logcounts[i,-zero])
  test <- cor.test(distance_hH[-zero], hH_logcounts[i,-zero])
  p_value_hH[i] <- test$p.value
}

names(cor_hH) <- rownames(rowData(spe_hH))
names(p_value_hH) <- rownames(rowData(spe_hH))

hist(cor_hH)
hist(abs(cor_hH))

# geni correlati (p.value significativo)
padjust_bf <- p.adjust(p_value_hH, method="fdr")
p_sig_bf <- which(padjust_bf < 0.01)
pv_sig_bf <- padjust_bf[p_sig_bf]

# AFTER batch effect removal
BC_hH_logcounts <- assays(spe_hH_noNA)$BClogcounts
BC_hH_logcounts <- as.matrix(BC_hH_logcounts)
distance_noNA <- spe_hH_noNA$inj_site_distance
cor_hH_noNA <- c()
p_value_hH_noNA <- c()

for(i in 1:nrow(BC_hH_logcounts)) {
  zero <- which(BC_hH_logcounts[i,] == 0)
  if(!isEmpty(zero)){
    cor_hH_noNA[i] <- cor(distance_noNA[-zero], BC_hH_logcounts[i,-zero])
    test <- cor.test(distance_noNA[-zero], BC_hH_logcounts[i,-zero])
    p_value_hH_noNA[i] <- test$p.value
  }else{
    if(isEmpty(zero)){
      cor_hH_noNA[i] <- cor(distance_noNA, BC_hH_logcounts[i,])
      test <- cor.test(distance_noNA, BC_hH_logcounts[i,])
      p_value_hH_noNA[i] <- test$p.value
    }
  }
}  

names(cor_hH_noNA) <- rownames(rowData(spe_hH_noNA))
names(p_value_hH_noNA) <- rownames(rowData(spe_hH_noNA))

hist(cor_hH_noNA)
hist(abs(cor_hH_noNA))

# geni correlati (p.value significativo)
padjust_af <- p.adjust(p_value_hH_noNA, method="fdr")
p_sig_af <- which(padjust_af < 0.01)
pv_sig_af <- padjust_af[p_sig_af]

write.csv(pv_sig_af, "Analysis - heme Hpx/geni_pvsig_hemeHpx.csv")
write.csv(cor_hH_noNA, "Analysis - heme Hpx/correlazioni_hemeHpx.csv")
