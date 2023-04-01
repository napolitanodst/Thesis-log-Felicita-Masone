# Extract significant genes for all modeling results

system.time(
  sig_genes <-
    sig_genes_extract_all(
      n = nrow(sce_layer),
      modeling_results = modeling_results,
      sce_layer = sce_layer
    )
)

# Extract the top n significant genes from the layer level modeling results

system.time(
  top5_sig_genes <- 
    sig_genes_extract(
      n = 5,
      modeling_results = modeling_results,
      model_type = names(modeling_results)[1],
      reverse = FALSE,
      sce_layer = sce_layer
    )
)
