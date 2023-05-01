# Extract significant genes

# Before BC
pval_bf <- rowData(spe)$pvalue_before
cor_bf <- rowData(spe)$cor_before

upg_bf <- which(cor_bf >= 0.35 & pval_bf < 0.01)
dng_bf <- which(cor_bf <= -0.35 & pval_bf < 0.01)

col <- rep("gray", nrow(genes))
col[upg_bf] <- "lightsalmon"
col[dng_bf] <- "seagreen1"

gene_name <- rowData(spe)$gene_name
nom <- c(gene_name[dng_bf], gene_name[upg_bf])

plot(x = cor_bf,
     y = -log(pval_bf),
     xlab = "Correlation",
     ylab = "-log p.values",
     type = "p",
     col = col,
     pch = 20,
     main = "Before batch correction"
)
abline(v=-0.35, col="black", lty=2)
abline(v=0.35, col="black", lty=2)
abline(h=-log(0.01), col="black", lty=2)
x = c(cor_bf[dng_bf], cor_bf[upg_bf])
y = c(pval_bf[dng_bf], pval_bf[upg_bf])
text(x, -log(y), 
     labels = nom,
     cex = 0.4)


# After BC
pval_af <- rowData(spe)$pvalue_after
cor_af <- rowData(spe)$cor_after

upg_af <- which(cor_af >= 0.35 & pval_af <= 0.01)
dng_af <- which(cor_af <= -0.35 & pval_af <= 0.01)

col2 <- rep("gray", nrow(genes))
col2[upg_af] <- "lightsalmon"
col2[dng_af] <- "seagreen1"

# vettore nomi
nom2 <- c(gene_name[dng_af], gene_name[upg_af])

plot(x = cor_af,
     y = -log(pval_af),
     xlab = "Correlation",
     ylab = "-log p.values",
     type = "p",
     col = col2,
     pch = 20,
     main = "After batch correction"
)
abline(v=-0.35, col="black", lty=2)
abline(v=0.35, col="black", lty=2)
abline(h=-log(0.01), col="black", lty=2)
x = c(cor_af[dng_af], cor_af[upg_af])
y = c(pval_af[dng_af], pval_af[upg_af])
text(x, -log(y), 
     labels = nom2,
     cex = 0.4)
