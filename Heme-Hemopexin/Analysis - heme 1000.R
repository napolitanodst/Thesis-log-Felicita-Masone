# SAMPLE: heme 1000 

# Import data as SpatialExperiment object ----
library(SpatialExperiment)
library(rtracklayer)
library(lobstr)

spe <- SpatialExperiment::read10xVisium(
  samples = "heme_1000",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# continuous variable containing total number of counts for each sample prior to filtering any genes
spe$sum_umi <- colSums(counts(spe))
# continuous variable containing the number of genes that have at least 1 count
spe$sum_gene <- colSums(counts(spe) > 0)

# gene annotation - reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/
gtf <-                                           
  rtracklayer::import(
    "gencode.vM23.annotation.gtf"
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

# Anatomy info and injection site coords ----
library(readr)

heme_1000_anatomy <- read_csv("heme_1000/heme_1000_anatomy.csv")
heme_1000_injection_site <- read_csv("heme_1000/heme_1000_injection_site.csv")

heme_1000_ana_inj <- merge(heme_1000_anatomy, heme_1000_injection_site, by = "Barcode")
remove(heme_1000_anatomy, heme_1000_injection_site)

spe$anatomy <- heme_1000_ana_inj$anatomy
spe$inj_site <- heme_1000_ana_inj$injection_site

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
  heme_1000_positions[i,7] <- sqrt((heme_1000_positions[i,"x"] - 52)^2 + (heme_1000_positions[i,"y"] - 70)^2)
}

colData(spe)$inj_site_distance <- heme_1000_positions$distance

# Quality Control and data filtering ----
library(scater)
library(ggspavis)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# histogram of numbers of expressed genes
hist(colData(spe)$detected, breaks = 20)

# Remove mitochondrial genes
spe <- spe[!is_mito, ]

# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(heme_1000_positions[,"anatomy"])))  
spe <- spe[,-NA_spot, drop = FALSE]

# remove not expressed genes
no_expr <- which(rowSums(counts(spe)) == 0)           
length(no_expr) / nrow(spe) * 100                     
spe <- spe[-no_expr, , drop = FALSE]

summary(spe$sum)

if (any(spe$sum == 0)) {
  spots_no_counts <- which(spe$sum == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

# Check the number of spots in which a gene is expressed and removal of genes expressed in less than 3 spots

h1000_counts <- assays(spe)$counts
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
spe <- spe[-no_rel, , drop = FALSE]
heme_1000_positions <- heme_1000_positions[-NA_spot,]

# Normalization
spe <- scuttle::logNormCounts(spe)

# Batch Correction on anatomy cluster ---- 
library(sva)

h1000_counts <- assays(spe)$counts
h1000_counts <- as.matrix(h1000_counts)

batch <- c()
for(i in 1:nrow(heme_1000_positions)){
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
heme_1000_positions$batch <- batch

## The combat_seq function also take "factors" as input, so another way to correct
## for anatomic batch is the following:
#  > batch <- as.factor(spe$anatomy)

BCcounts <- ComBat_seq(counts = h1000_counts,
                              batch = batch,
                              group = NULL,
                              covar_mod = NULL,
                              full_mod = FALSE)

assays(spe)$BCcounts <- BCcounts

spe <- scuttle::logNormCounts(x = spe,
                              assay.type = "BCcounts",
                              name = "BClogcounts")

# Reduced dimension: UMAP ----
library(scater)

# No corrected data UMAP
spe <- runUMAP(spe, exprs_values = "logcounts", name = "UMAP_noBC")
plotReducedDim(object = spe, dimred = "UMAP_noBC", colour_by = "anatomy")

# Batch Corrected data UMAP
spe <- runUMAP(spe, exprs_values = "BCcounts", name = "UMAP_BC")
plotReducedDim(object = spe, dimred = "UMAP_BC", colour_by = "anatomy")

spe <- runUMAP(spe, exprs_values = "BClogcounts", name = "UMAP_BClog")
plotReducedDim(object = spe, dimred = "UMAP_BClog", colour_by = "anatomy")

# Correlation and p.values ----

# BEFORE batch effect removal
counts_log <- assays(spe)$logcounts
counts_log <- as.matrix(counts_log)
distance_h <- spe$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_b <- c()
pval_b <- c()

for(i in 1:nrow(counts_log)) {
  zero <- which(counts_log[i,] == 0)
  if(!isEmpty(zero)){
    cor_b[i] <- cor(distance_h[-zero], counts_log[i,-zero], method = "spearman")
    test <- cor.test(distance_h[-zero], counts_log[i,-zero], method = "spearman")
    pval_b[i] <- test$p.value
  }else{
    if(isEmpty(zero)){
      cor_b[i] <- cor(distance_h, counts_log[i,], method = "spearman")
      test <- cor.test(distance_h, counts_log[i,], method = "spearman")
      pval_b[i] <- test$p.value
    }
  }
}

names(cor_b) <- rownames(rowData(spe))
names(pval_b) <- rownames(rowData(spe))

hist(cor_b)
hist(abs(cor_b))

# Correzione pvalue
padjust_bf <- p.adjust(pval_b, method="fdr")

# AFTER batch effect removal
BC_log <- assays(spe)$BClogcounts
BC_log <- as.matrix(BC_log)
cor_a <- c()
pval_a <- c()

for(i in 1:nrow(BC_log)) {
  zero <- which(BC_log[i,] == 0)
  if(!isEmpty(zero)){
    cor_a[i] <- cor(distance_h[-zero], BC_log[i,-zero], method = "spearman")
    test <- cor.test(distance_h[-zero], BC_log[i,-zero], method = "spearman")
    pval_a[i] <- test$p.value
  }else{
    if(isEmpty(zero)){
      cor_a[i] <- cor(distance_h, BC_log[i,], method = "spearman")
      test <- cor.test(distance_h, BC_log[i,], method = "spearman")
      pval_a[i] <- test$p.value
    }
  }
}

names(cor_a) <- rownames(rowData(spe))
names(pval_a) <- rownames(rowData(spe))

hist(cor_a)
hist(abs(cor_a))

# correzione pvalue
padjust_af <- p.adjust(pval_a, method="fdr")

# Add informations to spe
rowData(spe)$cor_before <- cor_b
rowData(spe)$pvalue_before <- padjust_bf

rowData(spe)$cor_after <- cor_a
rowData(spe)$pvalue_after <- padjust_af

# Compute teh number of 0s values
h1000_counts <- assays(spe)$counts
h1000_counts <- as.matrix(h1000_counts)
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

rowData(spe)$zero_x_gene <- zero_x_gene
rowData(spe)$no_zero_x_gene <- no_zero_x_gene

# Rank p.values
rank_bf <- rank(padjust_bf, ties.method = "min")
rank_af <- rank(padjust_af,  ties.method = "min")
rowData(spe)$rank_before <- rank_bf
rowData(spe)$rank_after <- rank_af
