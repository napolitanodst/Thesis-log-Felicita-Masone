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

# Quality Control and data filtering ----
library(scater)
library(ggspavis)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_h1000)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe_h1000 <- addPerCellQC(spe_h1000, subsets = list(mito = is_mito))

# histogram of numbers of expressed genes
hist(colData(spe_h1000)$detected, breaks = 20)

# Remove mitochondrial genes
spe_h1000 <- spe_h1000[!is_mito, ]

# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(heme_1000_positions[,5])))  
spe_h1000 <- spe_h1000[,-NA_spot, drop = FALSE]

# remove not expressed genes
no_expr <- which(rowSums(counts(spe_h1000)) == 0)           
length(no_expr) / nrow(spe_h1000) * 100                     
spe_h1000 <- spe_h1000[-no_expr, , drop = FALSE]

summary(spe_h1000$sum)

if (any(spe_h1000$sum == 0)) {
  spots_no_counts <- which(spe_h1000$sum == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe_h1000) * 100)
  spe_h1000 <- spe_h1000[, -spots_no_counts, drop = FALSE]
}

# Check the number of spots in which a gene is expressed and removal of genes expressed in less than 3 spots

h1000_counts <- assays(spe_h1000)$counts
h1000_counts <- as.matrix(h1000_counts)

no_rel <- c()
n <- 0

for(i in 1:nrow(h1000_counts)){
  for(j in 1:ncol(h1000_counts)){
    if(h1000_counts[i,j] != 0){
      n = n + 1
    }
  }
  if(n <= 2){
    no_rel = c(no_rel, i)
  }
  n = 0
}

remove(h1000_counts)
spe_h1000 <- spe_h1000[-no_rel, , drop = FALSE]
heme_1000_positions <- heme_1000_positions[-NA_spot,]

# Normalization
spe_h1000 <- scuttle::logNormCounts(spe_h1000)

# Batch Correction on anatomy cluster ---- 
library(sva)

h1000_counts <- assays(spe_h1000)$counts
h1000_counts <- as.matrix(h1000_counts)

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

BCcounts <- ComBat_seq(counts = h1000_counts,
                              batch = batch,
                              group = NULL,
                              covar_mod = NULL,
                              full_mod = FALSE)

assays(spe_h1000)$BC <- BCcounts

spe_h1000 <- logNormCounts(x = spe_h1000,
                            assay.type = "BC",
                            name = "BClogcounts")

# Reduced dimension: UMAP ----
library(scater)

# No corrected data UMAP
spe_h1000 <- runUMAP(spe_h1000, exprs_values = "logcounts", name = "UMAP_noBC")
plotReducedDim(spe_h1000, dimred = "UMAP_noBC", colour_by = "anatomy")

# Batch Corrected data UMAP
spe_h1000 <- runUMAP(spe_h1000, exprs_values = "BC", name = "UMAP_BC")
plotReducedDim(spe_h1000_noNA, dimred = "UMAP_BC", colour_by = "anatomy")

# Correlation and p.values ----

# BEFORE batch effect removal
h1000_logcounts <- assays(spe_h1000)$logcounts
h1000_logcounts <- as.matrix(h1000_logcounts)
distance_h <- spe_h1000$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_h <- c()
p_value_h <- c()

for(i in 1:nrow(h1000_logcounts)) {
  zero <- which(h1000_logcounts[i,] == 0)
  if(!isEmpty(zero)){
    cor_h[i] <- cor(distance_h[-zero], h1000_logcounts[i,-zero])
    test <- cor.test(distance_h[-zero], h1000_logcounts[i,-zero])
    p_value_h[i] <- test$p.value
  }else{
    if(isEmpty(zero)){
      cor_h[i] <- cor(distance_h, h1000_logcounts[i,])
      test <- cor.test(distance_h, h1000_logcounts[i,])
      p_value_h[i] <- test$p.value
    }
  }
}

names(cor_h) <- rownames(rowData(spe_h1000))
names(p_value_h) <- rownames(rowData(spe_h1000))

hist(cor_h)
hist(abs(cor_h))

# correlated genes before batch correction
padjust_bf <- p.adjust(p_value_h, method="fdr")
pv_sig_bf <- padjust_bf[which(padjust_bf < 0.01)]

# AFTER batch effect removal
BC_h_logcounts <- assays(spe_h1000)$BClogcounts
BC_h_logcounts <- as.matrix(BC_h_logcounts)
cor_h_BC <- c()
p_value_h_BC <- c()

for(i in 1:nrow(BC_h_logcounts)) {
  zero <- which(BC_h_logcounts[i,] == 0)
  if(!isEmpty(zero)){
    cor_h_BC[i] <- cor(distance_h[-zero], BC_h_logcounts[i,-zero])
    test <- cor.test(distance_h[-zero], BC_h_logcounts[i,-zero])
    p_value_h_BC[i] <- test$p.value
  }else{
    if(isEmpty(zero)){
      cor_h_BC[i] <- cor(distance_h, BC_h_logcounts[i,])
      test <- cor.test(distance_h, BC_h_logcounts[i,])
      p_value_h_BC[i] <- test$p.value
     }
   }
 }

names(cor_h_BC) <- rownames(rowData(spe_h1000))
names(p_value_h_BC) <- rownames(rowData(spe_h1000))

hist(cor_h_BC)
hist(abs(cor_h_BC))

# correlated genes after batch correction
padjust_af <- p.adjust(p_value_h_BC, method="fdr")
pv_sig_af <- padjust_af[which(padjust_af < 0.01)]

cor_pv_sig_h <- cor_h_BC[which(padjust_af < 0.01)]

# Add all these informations to spe ----
rowData(spe_h1000)$cor_before <- cor_h
rowData(spe_h1000)$pvalue_before <- padjust_bf
rowData(spe_h1000)$cor_after <- cor_h_BC
rowData(spe_h1000)$pvalue_after <- padjust_af

# Compute teh number of 0s values
zero_x_gene <- c()
n <- 0

for(i in 1:nrow(h1000_counts)){
  for(j in 1:ncol(h1000_counts)){
    if(h1000_counts[i,j] == 0){
      n = n + 1
    }
  }
  zero_x_gene = c(zero_x_gene, n)
  n = 0
}

no_zero_x_gene <- 2310 - zero_x_gene 

rowData(spe_h1000)$zero_x_gene <- zero_x_gene
rowData(spe_h1000)$no_zero_x_gene <- no_zero_x_gene

# Rank p.values
rank_bf <- rank(padjust_bf)
rank_af <- rank(padjust_af)
rowData(spe_h1000)$rank_before <- rank_bf
rowData(spe_h1000)$rank_after <- rank_af
