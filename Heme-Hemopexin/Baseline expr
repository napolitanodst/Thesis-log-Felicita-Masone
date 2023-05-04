## Baseline

non0s_x_gene <- rowData(spe_h1000)$no_zero_x_gene

# before
counts_matrix <- assays(spe)$counts
counts_matrix <- as.matrix(counts_matrix)
average_expr_before <- rowSums(counts_matrix)/non0s_x_gene
mean_expr_before <- mean(average_expr_before)

dist_mean_before <- c()
for(i in 1:length(average_expr_before)){
  dist_mean_before[i] = abs(average_expr_before[i] - mean_expr_before)
}

rank_baseline_before <- rank(dist_mean_before, ties.method = "min")

rowData(spe)$baseline_before <- average_expr_before
rowData(spe)$rank_baseline_before <- rank_baseline_before

# after
BCcounts <- assays(spe)$BCcounts
BCcounts <- as.matrix(BCcounts)
average_expr_after <- rowSums(BCcounts)/non0s_x_gene
mean_expr_after <- mean(average_expr_after)

dist_mean_after <- c()
for(i in 1:length(average_expr_after)){
  dist_mean_after[i] = abs(average_expr_after[i] - mean_expr_after)
}

rank_baseline_after <- rank(dist_mean_after, ties.method = "min")

rowData(spe)$baseline_after <- average_expr_after
rowData(spe)$rank_baseline_after <- rank_baseline_after
