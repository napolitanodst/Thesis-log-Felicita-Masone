# txt file containing the genes ID stored in rowDAta(spe)
write.table(rowData(spe)$gene_search, "C:/Users/mason/OneDrive/Documenti/genes_ID.txt" )

# Sample Spatial Gene visualization
vis_gene(
  spe = spe,
  sampleid = "151507",
  geneid = "MED11; ENSG00000161920",  ## A character(1) specifyng the gene ID stored in rowData(spe)$gene_search or a continuous variable stored in colData(spe) to visualize. Also rownames(spe) can be used to search for the gene ID
  spatial = TRUE,
  assayname = "logcounts",
  image_id = "lowres",
  alpha = 1,
  point_size = 2,
  ... = " "
)


