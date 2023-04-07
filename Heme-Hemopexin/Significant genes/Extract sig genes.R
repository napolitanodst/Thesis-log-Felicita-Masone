# Extract significant genes

# Before BC
pvalue_bf <- rowData(spe_h1000)$pvalue_before
cor_bf <- rowData(spe_h1000)$cor_before

upg_bf <- which(cor_bf >= 0.35 & pvalue_bf < 0.01)
dng_bf <- which(cor_bf <= -0.35 & pvalue_bf < 0.01)

col <- rep("gray", nrow(tab))
col[upg_bf] <- "pink"
col[dng_bf] <- "green"

gene_name <- rowData(spe_h1000)$gene_name
nom <- c(gene_name[dng_bf], gene_name[upg_bf])

plot(x = cor_bf,
     y = -log(pvalue_bf),
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
y = c(pvalue_bf[dng_bf], pvalue_bf[upg_bf])
text(x, -log(y), 
     labels = nom,
     cex = 0.4)


# After BC
pvalue_af <- rowData(spe_h1000)$pvalue_after
cor_af <- rowData(spe_h1000)$cor_after

upg_af <- which(cor_af >= 0.35 & pvalue_af <= 0.01)
dng_af <- which(cor_af <= -0.35 & pvalue_af <= 0.01)

col2 <- rep("gray", nrow(tab))
col2[upg_af] <- "pink"
col2[dng_af] <- "green"

# vettore nomi
nom2 <- c(gene_name[dng_af], gene_name[upg_af])

plot(x = cor_af,
     y = -log(pvalue_af),
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
y = c(pvalue_af[dng_af], pvalue_af[upg_af])
text(x, -log(y), 
     labels = nom2,
     cex = 0.4)
