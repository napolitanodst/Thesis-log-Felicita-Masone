## SVG 
library(spatialDE)

# BEFORE BATCH CORRECTION
# preparing data for input
x_cord <- spe$array_row
y_cord <- spe$array_col
spots <- colnames(spe)
sample_info <- matrix(data = NA, nrow = 2310, ncol = 2)
colnames(sample_info) <- c("x", "y")
rownames(sample_info) <- spots
sample_info[,1] <- x_cord
sample_info[,2] <- y_cord
sample_info <- as.data.frame(sample_info)
gene_id <- rowData(spe)$gene_id
expr_counts <- assays(spe)$counts
expr_counts <- as.matrix(expr_counts)
total_counts <- colSums(expr_counts)
sample_info$total_counts <- total_counts
coordinates <- sample_info[,-3]

# spatialDE pipeline
expr_norm <- stabilize(expr_counts)
expr_resid <- regress_out(expr_norm, sample_info = sample_info)
results <- spatialDE::run(expr_resid, coordinates = coordinates)

de_results <- results[results$qval < 0.01, ] 

ms_results <- model_search(
  expr_resid,
  coordinates = coordinates,
  de_results = de_results
)

SVG_before <- results[,c("g", "LLR", "pval")]

# Pvalues are computed on the LLR statistic, so we rank these values
# cor_pLLR <- cor(SVG_before$LLR, SVG_before$pval)
# > -0.2961986
# spearman -0.9997527
rank_llr <- rank(SVG_before$LLR, ties.method = "min")

# Ranks performed on LLR classify genes from least significant to most significant so it is necessary to invert ranks in order to perform method comparison analyses
SVG_before$rank_llr <- 16187 - rank_llr
rownames(SVG_before) <- c(1:16186)

# AFTER BATCH CORRECTION
sample_info_BC <- sample_info[,-3]
BCcounts <- assays(spe)$BCcounts
BCcounts <- as.matrix(BCcounts)
total_counts_BC <- colSums(BCcounts)
sample_info_BC$total_counts <- total_counts_BC

expr_norm_BC <- stabilize(BCcounts)
expr_resid_BC <- regress_out(expr_norm_BC, sample_info = sample_info_BC)
results_BC <- spatialDE::run(expr_resid_BC, coordinates = coordinates)

de_results_BC <- results_BC[results_BC$qval < 0.01, ] 

ms_results_BC <- model_search(
  expr_resid_BC,
  coordinates = coordinates,
  de_results = de_results_BC
)

SVG_BC <- results_BC[,c("g", "LLR", "pval")]
rownames(SVG_BC) <- c(1:16186)

rank_llr_BC <- rank(SVG_BC$LLR, ties.method = "min")
SVG_BC$rank_llr <- 16187 - rank_llr_BC
rownames(SVG_BC) <- c(1:16186)

SVG <- merge(SVG_before, SVG_BC, by = "gene_id")

# ground truth gene positions in SVG
SVG_up <- c()
up_id <- heme_upgenes$id
for (i in 1:62) {
  SVG_up[i] <- which(SVG$g == up_id[i])  
}

SVG_dn <- c()
dn_id <- heme_dngenes$id
for (i in 1:6) {
  SVG_dn[i] <- which(SVG$g == dn_id[i])  
}
