library(spatialLIBD)

## Download the Human DLPFC Visium data from LIBD

ehub <- ExperimentHub::ExperimentHub()

spe <- fetch_data(type = "spe", eh = ehub)
sce <- fetch_data(type = "sce", eh = ehub)
sce_layer <- fetch_data(type = "sce_layer", eh = ehub)
modeling_results <- fetch_data("modeling_results", eh = ehub)