# SAMPLE: control 

# Import data as SpatialExperiment object ----
library(SpatialExperiment)
library(rtracklayer)
library(lobstr)

spe_c <- SpatialExperiment::read10xVisium(
  samples = "control",
  sample_id = "control",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# continuous variable containing total number of counts for each sample prior to filtering any genes
spe_c$sum_umi <- colSums(counts(spe_c))
# continuous variable containing the number of genes that have at least 1 count
spe_c$sum_gene <- colSums(counts(spe_c) > 0)

spe_c <- scuttle::logNormCounts(spe_c)

# gene annotation - reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/
gtf <-                                           
  rtracklayer::import(
    "gencode.vM23.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    
names(gtf) <- gtf$gene_id                        

match_genes <- match(rownames(spe_c), gtf$gene_id) 
table(is.na(match_genes))
spe_c <- spe_c[!is.na(match_genes), ]                
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  
rowRanges(spe_c) <- gtf[match_genes]               

rowData(spe_c)$gene_search <- paste0(
  rowData(spe_c)$gene_name, "; ", rowData(spe_c)$gene_id
)

# Filter spe from mitochondrial and not expressed genes
is_mito <- which(seqnames(spe_c) == "chrM")
spe_c$expr_chrM <- colSums(counts(spe_c)[is_mito, , drop = FALSE])
spe_c$expr_chrM_ratio <- spe_c$expr_chrM / spe_c$sum_umi

no_expr <- which(rowSums(counts(spe_c)) == 0)           
length(no_expr) / nrow(spe_c) * 100                     
spe_c <- spe_c[-no_expr, , drop = FALSE]

summary(spe_c$sum_umi)

if (any(spe_c$sum_umi == 0)) {
  spots_no_counts <- which(spe_c$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe_c) * 100)
  spe_c <- spe_c[, -spots_no_counts, drop = FALSE]
}

# Anatomy info and injection site coords ----
library(readr)

control_anatomy <- read_csv("control/control_anatomy.csv")
control_injection_site <- read_csv("control/control_injection_site.csv")

control_ana_inj <- merge(control_anatomy, control_injection_site, by = "Barcode")
remove(control_anatomy, control_injection_site)

spe_c$anatomy <- control_ana_inj$anatomy
spe_c$inj_site <- control_ana_inj$injection_site

# Distances from injection site ----

tissue_positions_list_control <- read_csv("control/outs/spatial/tissue_positions_list.csv", col_names = FALSE)
names(tissue_positions_list_control) <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
tissue_positions_list_control <- subset(tissue_positions_list_control, in_tissue == 1)
tissue_positions_list_control <- tissue_positions_list_control[,-c(5,6)]
control_positions <- merge(tissue_positions_list_control, control_ana_inj, by = "Barcode")
remove(tissue_positions_list_control)

# calcolo la distanza utilizzando le coordinate dello spot centrale del sito di iniezione 
# x = 44 ; y = 70
control_positions$distance <- NA
for(i in 1:nrow(control_positions)){
  control_positions[i,7] <- sqrt((control_positions[i,3] - 44)^2 + (control_positions[i,4] - 70)^2)
}

colData(spe_c)$inj_site_distance <- control_positions$distance

# Batch Correction on anatomy cluster ---- 
library(sva)

expr_counts <- assays(spe_c)$counts
expr_counts <- as.matrix(expr_counts)

batch <- c()
for(i in 1:nrow(control_positions)){
  if(is.na(control_positions[i,5])){
    batch[i] <- NA
  }else{
    if(control_positions[i,5] == "caudate_putamen"){
      batch[i] <- 1
    }else{
      if(control_positions[i,5] == "cortex"){
        batch[i] <- 2
      }else{
        if(control_positions[i,5] == "thalamus"){
          batch[i] <- 3
        }else{
          if(control_positions[i,5] == "globus_pallidus"){
            batch[i] <- 4
          }else{
            if(control_positions[i,5] == "plexus"){
              batch[i] <- 5
            }else{
              if(control_positions[i,5] == "hypothalamus"){
                batch[i] <- 6
              }else{
                if(control_positions[i,5] == "corpus_callosum"){
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
control_positions$batch <- batch

# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(batch)))  
batch_noNA <- batch[!is.na(batch)]
expr_counts_noNA <- expr_counts[,-NA_spot] 

spe_c_noNA <- spe_c[,-NA_spot]

BCmatrix_noNA_c <- ComBat_seq(counts = expr_counts_noNA,
                            batch = batch_noNA,
                            group = NULL,
                            covar_mod = NULL,
                            full_mod = FALSE)

assays(spe_c_noNA)$BC <- BCmatrix_noNA_c

# Reduced dimension: UMAP ----
library(scater)

# No corrected data UMAP
spe_c <- runUMAP(spe_c, exprs_values = "logcounts")
reducedDimNames(spe_c)
head(reducedDim(spe_c))
plotUMAP(spe_c, colour_by = "anatomy")

# Batch Corrected data UMAP
spe_c_noNA <- runUMAP(spe_c_noNA, exprs_values = "BC")
reducedDimNames(spe_c_noNA)
head(reducedDim(spe_c_noNA))
plotUMAP(spe_c_noNA, colour_by = "anatomy")

# Correlation and p.values ----

# BEFORE batch effect removal
expr_log_c <- assays(spe_c)$logcounts
expr_log_c <- as.matrix(expr_log_c)
distance_c <- spe_c$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_c <- c()
p_value_c <- c()

for(i in 1:nrow(expr_log_c)) {
  cor_c[i] <- cor(distance_c, expr_log_c[i,])
  test <- cor.test(distance_c, expr_log_c[i,])
  p_value_c[i] <- test$p.value
}

names(cor_c) <- rownames(rowData(spe_c))
names(p_value_c) <- rownames(rowData(spe_c))

hist(cor_c)
hist(abs(cor_c))

# geni correlati (p.value significativo)
padjust_bf <- p.adjust(p_value_c, method="fdr")
p_sig_bf <- which(padjust_bf < 0.01)
pv_sig_bf <- padjust_bf[p_sig_bf]

# AFTER batch effect removal
distance_noNA <- spe_c_noNA$inj_site_distance
cor_c_noNA <- c()
p_value_c_noNA <- c()

for(i in 1:nrow(BCmatrix_noNA_c)) {
  cor_c_noNA[i] <- cor(distance_noNA, BCmatrix_noNA_c[i,])
  test <- cor.test(distance_noNA, BCmatrix_noNA_c[i,])
  p_value_c_noNA[i] <- test$p.value
}

names(cor_c_noNA) <- rownames(rowData(spe_c_noNA))
names(p_value_c_noNA) <- rownames(rowData(spe_c_noNA))

hist(cor_c_noNA)
hist(abs(cor_c_noNA))

# geni correlati (p.value significativo)
padjust_af <- p.adjust(p_value_c_noNA, method="fdr")
p_sig_af <- which(padjust_af < 0.01)
pv_sig_af <- padjust_af[p_sig_af]

write.csv(pv_sig_af, "Analysis - control/geni_pvsig.csv")
write.csv(cor_c_noNA, "Analysis - control/correlazioni.csv")
