# 21/03/2023 

# Quality Control
library(scater)
library(ggspavis)

# --- Control ----

# calculateQCMetrics : calculateQCMetrics is succeeded by perCellQCMetrics and perFeatureQCMetrics;
#                      normalize is succeeded by logNormCounts;

spe_c <- addPerCellQCMetrics(spe_c)

## Library size - sum_umi

# histogram of library sizes : Library size represents the total sum of UMI counts per spot
hist(colData(spe_c)$sum_umi, breaks = 20)

#We set a relatively arbitrary threshold of 600 UMI counts per spot, and then check the number of spots below this threshold.
# select QC threshold for library size
qc_lib_size_control <- colData(spe_c)$sum_umi < 600
table(qc_lib_size_control)
colData(spe_c)$qc_lib_size <- qc_lib_size_control

# check spatial pattern of discarded spots
plotQC(spe_c, type = "spots", 
       discard = "qc_lib_size")

## Number of expressed features - detected

# histogram of numbers of expressed genes 
hist(colData(spe_c)$detected, breaks = 20)

# select QC threshold for number of expressed genes
qc_detected_control <- colData(spe_c)$detected < 400
table(qc_detected_control)
colData(spe_c)$qc_detected <- qc_detected_control

# check spatial pattern of discarded spots
plotQC(spe_c, type = "spots", 
       discard = "qc_detected")


# --- Heme 1000 ----

# calculate QC metrics
spe_h1000 <- addPerCellQCMetrics(spe_h1000)

## Library size - sum_umi

# histogram of library sizes : Library size represents the total sum of UMI counts per spot
hist(colData(spe_h1000)$sum_umi, breaks = 20)

#We set a relatively arbitrary threshold of 600 UMI counts per spot, and then check the number of spots below this threshold.
# select QC threshold for library size
qc_lib_size_heme1000 <- colData(spe_h1000)$sum_umi < 600
table(qc_lib_size_heme1000)
colData(spe_h1000)$qc_lib_size <- qc_lib_size_heme1000

# check spatial pattern of discarded spots
plotQC(spe_h1000, type = "spots", 
       discard = "qc_lib_size")

## Number of expressed features - detected

# histogram of numbers of expressed genes 
hist(colData(spe_h1000)$detected, breaks = 20)

# select QC threshold for number of expressed genes
qc_detected_heme1000 <- colData(spe_h1000)$detected < 400
table(qc_detected_heme1000)
colData(spe_h1000)$qc_detected <- qc_detected_heme1000

# check spatial pattern of discarded spots
plotQC(spe_h1000, type = "spots", 
       discard = "qc_detected")
