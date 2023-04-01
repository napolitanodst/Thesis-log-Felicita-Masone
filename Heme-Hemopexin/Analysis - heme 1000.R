# SAMPLE: heme 1000 

# Import data as SpatialExperiment object ----
library(SpatialExperiment)
library(rtracklayer)
library(lobstr)

spe_h1000 <- SpatialExperiment::read10xVisium(
  samples = "heme_1000",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# continuous variable containing total number of counts for each sample prior to filtering any genes
spe_h1000$sum_umi <- colSums(counts(spe_h1000))
# continuous variable containing the number of genes that have at least 1 count
spe_h1000$sum_gene <- colSums(counts(spe_h1000) > 0)

spe_h1000 <- scuttle::logNormCounts(spe_h1000)

# gene annotation - reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/
gtf <-                                           
  rtracklayer::import(
    "gencode.vM23.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    
names(gtf) <- gtf$gene_id                        

match_genes <- match(rownames(spe_h1000), gtf$gene_id) 
table(is.na(match_genes))
spe_h1000 <- spe_h1000[!is.na(match_genes), ]                
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  
rowRanges(spe_h1000) <- gtf[match_genes]               

rowData(spe_h1000)$gene_search <- paste0(
  rowData(spe_h1000)$gene_name, "; ", rowData(spe_h1000)$gene_id
)

# Filter spe from mitochondrial and not expressed genes
is_mito <- which(seqnames(spe_h1000) == "chrM")
spe_h1000$expr_chrM <- colSums(counts(spe_h1000)[is_mito, , drop = FALSE])
spe_h1000$expr_chrM_ratio <- spe_h1000$expr_chrM / spe_h1000$sum_umi

no_expr <- which(rowSums(counts(spe_h1000)) == 0)           
length(no_expr) / nrow(spe_h1000) * 100                     
spe_h1000 <- spe_h1000[-no_expr, , drop = FALSE]

summary(spe_h1000$sum_umi)

if (any(spe_h1000$sum_umi == 0)) {
  spots_no_counts <- which(spe_h1000$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe_h1000) * 100)
  spe_h1000 <- spe_h1000[, -spots_no_counts, drop = FALSE]
}

# Anatomy info and injection site coords ----
library(readr)

heme_1000_anatomy <- read_csv("heme_1000/heme_1000_anatomy.csv")
heme_1000_injection_site <- read_csv("heme_1000/heme_1000_injection_site.csv")

heme_1000_ana_inj <- merge(heme_1000_anatomy, heme_1000_injection_site, by = "Barcode")
remove(heme_1000_anatomy, heme_1000_injection_site)

spe_h1000$anatomy <- heme_1000_ana_inj$anatomy
spe_h1000$inj_site <- heme_1000_ana_inj$injection_site

# Distances from injection site ----

tissue_positions_list_h1000 <- read_csv("heme_1000/outs/spatial/tissue_positions_list.csv", col_names = FALSE)
names(tissue_positions_list_h1000) <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
tissue_positions_list_h1000 <- subset(tissue_positions_list_h1000, in_tissue == 1)
tissue_positions_list_h1000 <- tissue_positions_list_h1000[,-c(5,6)]
heme_1000_positions <- merge(tissue_positions_list_h1000, heme_1000_ana_inj, by = "Barcode")
remove(tissue_positions_list_h1000)

# calcolo la distanza utilizzando le coordinate dello spot centrale del sito di iniezione 
# x = 52 ; y = 70
heme_1000_positions$distance <- NA
for(i in 1:nrow(heme_1000_positions)){
  heme_1000_positions[i,7] <- sqrt((heme_1000_positions[i,3] - 52)^2 + (heme_1000_positions[i,4] - 70)^2)
}

colData(spe_h1000)$inj_site_distance <- heme_1000_positions$distance

# Batch Correction on anatomy cluster ---- 
library(sva)

expr_counts <- assays(spe_h1000)$counts
expr_counts <- as.matrix(expr_counts)

batch <- c()
for(i in 1:nrow(heme_1000_positions)){
  if(is.na(heme_1000_positions[i,5])){
    batch[i] <- NA
  }else{
    if(heme_1000_positions[i,5] == "caudate_putamen"){
      batch[i] <- 1
    }else{
      if(heme_1000_positions[i,5] == "cortex"){
        batch[i] <- 2
      }else{
        if(heme_1000_positions[i,5] == "thalamus"){
          batch[i] <- 3
        }else{
          if(heme_1000_positions[i,5] == "globus_pallidus"){
            batch[i] <- 4
          }else{
            if(heme_1000_positions[i,5] == "plexus"){
              batch[i] <- 5
            }else{
              if(heme_1000_positions[i,5] == "hypothalamus"){
                batch[i] <- 6
              }else{
                if(heme_1000_positions[i,5] == "corpus_callosum"){
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
heme_1000_positions$batch <- batch

# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(batch)))  
batch_noNA <- batch[!is.na(batch)]
expr_counts_noNA <- expr_counts[,-NA_spot] 

spe_h1000_noNA <- spe_h1000[,-NA_spot]
spe_h1000_noNA$sum_umi <- colSums(counts(spe_h1000_noNA))
spe_h1000_noNA$sum_gene <- colSums(counts(spe_h1000_noNA) > 0)

BCmatrix_noNA_h <- ComBat_seq(counts = expr_counts_noNA,
                              batch = batch_noNA,
                              group = NULL,
                              covar_mod = NULL,
                              full_mod = FALSE)

assays(spe_h1000_noNA)$BC <- BCmatrix_noNA_h

# Reduced dimension: UMAP ----
library(scater)

# No corrected data UMAP
spe_h1000 <- runUMAP(spe_h1000, exprs_values = "logcounts")
reducedDimNames(spe_h1000)
head(reducedDim(spe_h1000))
plotUMAP(spe_h1000, colour_by = "anatomy")

# Batch Corrected data UMAP
spe_h1000_noNA <- runUMAP(spe_h1000_noNA, exprs_values = "BC")
reducedDimNames(spe_h1000_noNA)
head(reducedDim(spe_h1000_noNA))
plotUMAP(spe_h1000_noNA, colour_by = "anatomy")

# Correlation and p.values ----

# BEFORE batch effect removal
expr_log_h <- assays(spe_h1000)$logcounts
expr_log_h <- as.matrix(expr_log_h)
distance_h <- spe_h1000$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_h <- c()
p_value_h <- c()

for(i in 1:nrow(expr_log_h)) {
  cor_h[i] <- cor(distance_h, expr_log_h[i,])
  test <- cor.test(distance_h, expr_log_h[i,])
  p_value_h[i] <- test$p.value
}

names(cor_h) <- rownames(rowData(spe_h1000))
names(p_value_h) <- rownames(rowData(spe_h1000))

hist(cor_h)
hist(abs(cor_h))

# geni correlati (p.value significativo)
padjust_bf <- p.adjust(p_value_h, method="fdr")
p_sig_bf <- which(padjust_bf < 0.01)
pv_sig_bf <- padjust_bf[p_sig_bf]

# AFTER batch effect removal
distance_noNA <- spe_h1000_noNA$inj_site_distance
cor_h_noNA <- c()
p_value_h_noNA <- c()

for(i in 1:nrow(BCmatrix_noNA_h)) {
  cor_h_noNA[i] <- cor(distance_noNA, BCmatrix_noNA_h[i,])
  test <- cor.test(distance_noNA, BCmatrix_noNA_h[i,])
  p_value_h_noNA[i] <- test$p.value
}

names(cor_h_noNA) <- rownames(rowData(spe_h1000_noNA))
names(p_value_h_noNA) <- rownames(rowData(spe_h1000_noNA))

hist(cor_h_noNA)
hist(abs(cor_h_noNA))

# geni correlati (p.value significativo)
padjust_af <- p.adjust(p_value_h_noNA, method="fdr")
p_sig_af <- which(padjust_af < 0.01)
pv_sig_af <- padjust_af[p_sig_af]

write.csv(pv_sig_af, "Analysis - heme 1000/geni_pvsig.csv")
write.csv(cor_h_noNA, "Analysis - heme 1000/correlazioni.csv")
