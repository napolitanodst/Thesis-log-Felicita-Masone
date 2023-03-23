# --- Correlazione tra la distanza e l'espressione genica (logcounts) per
# ogni gene e relativi p.value

## Control ----

# Estraggo i valori di espressione dall'oggetto spe_c. La matrice dei conteggi è
# di tipo sparse matrix of class DelayedMatrix, perciò la rendo matrice di tipo matrix. 
# Nello slot assays oltre alla matrice dei conteggi c'è quella dei logcounts

expr_log_control <- assays(spe_c)$logcounts
expr_log_control <- as.matrix(expr_log_control)
control_distance <- spe_c$inj_site_distance

# Calcolo la correlazione tra la distanza e l'espressione per ogni gene e il relativo p.value
cor_control <- c()
p_value_control <- c()

i <- 1
while(i <= nrow(expr_log_control)) {
  cor_control[i] <- cor(control_distance, expr_log_control[i,])
  test <- cor.test(control_distance, expr_log_control[i,])
  p_value_control[i] <- test$p.value
  i <- i + 1
  if(i > nrow(expr_log_control)){rm(test)}
}

names(cor_control) <- rownames(rowData(spe_c))
names(p_value_control) <- rownames(rowData(spe_c))

## Heme 1000 ----

expr_log_heme1000 <- assays(spe_h1000)$logcounts
expr_log_heme1000 <- as.matrix(expr_log_heme1000)
heme1000_distance <- spe_h1000$inj_site_distance

cor_heme1000 <- c()
p_value_heme1000 <- c()

i <- 1
while(i <= nrow(expr_log_heme1000)) {
  cor_heme1000[i] <- cor(heme1000_distance, expr_log_heme1000[i,])
  test <- cor.test(heme1000_distance, expr_log_heme1000[i,])
  p_value_heme1000[i] <- test$p.value
  i <- i + 1
  if(i > nrow(expr_log_heme1000)){rm(test)}
}

names(cor_heme1000) <- rownames(rowData(spe_h1000))
names(p_value_heme1000) <- rownames(rowData(spe_h1000))
