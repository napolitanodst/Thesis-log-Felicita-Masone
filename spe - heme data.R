# Create a spe object with the Spatial transcriptome data from coronal mouse brain 
# sections after striatal injection of heme and heme-hemopexin. 

# Control: mouse brain without any treatment
spe <- SpatialExperiment::read10xVisium(
  samples = "control",
  sample_id = "control",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# reference genome GRCm38.p6 (mm10) https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz




