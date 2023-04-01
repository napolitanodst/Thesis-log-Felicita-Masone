# Sample Spatial Cluster visualization

vis_clus(
  spe = spe,
  sampleid = "151507",
  clustervar = "spatialLIBD",  ## A character(1) with the name of the colData(spe) column that has the cluster values
  colors = libd_layer_colors,
  spatial = TRUE,              ## A logical(1) indicating wheter to include the histology layer from geom_spatial()
  image_id = "lowres",         
  alpha = 1,                   ## A numeric(1) in the [0, 1] range that specifies the transparency level of the data on the spots
  point_size = 2,              ## A numeric(1) specifyng the size of the point
  ... = " "
)
