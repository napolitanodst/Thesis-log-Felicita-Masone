## Baseline

non0s_x_gene <- rowData(spe_h1000)$no_zero_x_gene

# before
counts_matrix <- assays(spe)$counts
counts_matrix <- as.matrix(counts_matrix)
baseline_before <- rowSums(counts_matrix)/non0s_x_gene
mean_expr_before <- mean(baseline_before)

dist_mean_before <- c()
for(i in 1:length(baseline_before)){
  dist_mean_before[i] = abs(baseline_before[i] - mean_expr_before)
}

rank_baseline_before <- rank(dist_mean_before, ties.method = "min")

rowData(spe)$baseline_before <- baseline_before
rowData(spe)$rank_baseline_before <- rank_baseline_before

# after
BCcounts <- assays(spe)$BCcounts
BCcounts <- as.matrix(BCcounts)
baseline_after <- rowSums(BCcounts)/non0s_x_gene
mean_expr_after <- mean(baseline_after)

dist_mean_after <- c()
for(i in 1:length(baseline_after)){
  dist_mean_after[i] = abs(baseline_after[i] - mean_expr_after)
}

rank_baseline_after <- rank(dist_mean_after, ties.method = "min")

rowData(spe)$baseline_after <- baseline_after
rowData(spe)$rank_baseline_after <- rank_baseline_after
