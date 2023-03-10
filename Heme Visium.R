## -- VISIUM LIEBER EXAMPLE with the Spatial transcriptome data from coronal mouse brain sections 
##    after striatal injection of heme and heme-hemopexin


library(tidyverse)
library(ggplot2)
library(Matrix)
library(Rmisc)
library(ggforce)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)


# --- Reading data ---

#1 Define samples
sample_names <- c("control", "heme_0030", "heme_0125", "heme_0500", "heme_1000", "hemeHpx_1000", "sham")

#2 Define your paths (Paths should be in the same order as the corresponding sample names)

image_paths <- c("input/input/control/spatial/tissue_lowres_image.png",
                 "input/input/heme_0030/spatial/tissue_lowres_image.png",
                 "input/input/heme_0125/spatial/tissue_lowres_image.png",
                 "input/input/heme_0500/spatial/tissue_lowres_image.png",
                 "input/input/heme_1000/spatial/tissue_lowres_image.png",
                 "input/input/hemeHpx_1000/spatial/tissue_lowres_image.png",
                 "input/input/sham/spatial/tissue_lowres_image.png"
)

scalefactor_paths <- c("input/input/control/spatial/scalefactors_json.json",
                       "input/input/heme_0030/spatial/scalefactors_json.json",
                       "input/input/heme_0125/spatial/scalefactors_json.json",
                       "input/input/heme_0500/spatial/scalefactors_json.json",
                       "input/input/heme_1000/spatial/scalefactors_json.json",
                       "input/input/hemeHpx_1000/spatial/scalefactors_json.json",
                       "input/input/sham/spatial/scalefactors_json.json"
)

tissue_paths <- c("input/input/control/spatial/tissue_positions_list.txt",
                  "input/input/heme_0030/spatial/tissue_positions_list.txt",
                  "input/input/heme_0125/spatial/tissue_positions_list.txt",
                  "input/input/heme_0500/spatial/tissue_positions_list.txt",
                  "input/input/heme_1000/spatial/tissue_positions_list.txt",
                  "input/input/hemeHpx_1000/spatial/tissue_positions_list.txt",
                  "input/input/sham/spatial/tissue_positions_list.txt"
)

anatomy_paths <- c("input/input/control/control_anatomy.csv",
                   "input/input/heme_0030/heme_0030_anatomy.csv",
                   "input/input/heme_0125/heme_0125_anatomy.csv",
                   "input/input/heme_0500/heme_0500_anatomy.csv",                    
                   "input/input/heme_1000/heme_1000_anatomy.csv",
                   "input/input/hemeHpx_1000/hemeHpx_1000_anatomy.csv",
                   "input/input/sham/sham_anatomy.csv"
)


matrix_paths <- c("input/input/control/filtered_feature_bc_matrix.h5",
                  "input/input/heme_0030/filtered_feature_bc_matrix.h5",
                  "input/input/heme_0125/filtered_feature_bc_matrix.h5",
                  "input/input/heme_0500/filtered_feature_bc_matrix.h5",                    
                  "input/input/heme_1000/filtered_feature_bc_matrix.h5",
                  "input/input/hemeHpx_1000/filtered_feature_bc_matrix.h5",
                  "input/input/sham/filtered_feature_bc_matrix.h5"
)

#3 Read in down sampled images (We need to determine the image height and width for proper plotting in the end)
images_cl <- list()

for (i in 1:length(sample_names)) {
  images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
  height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
  width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)

#4 Convert the images to grobs (This step provides compatibility with ggplot2)
grobs <- list()
for (i in 1:length(sample_names)) {
  grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width

images_tibble

#5 Read in Clusters
clusters <- list()
for (i in 1:length(sample_names)) {
  clusters[[i]] <- read.csv(anatomy_paths[i])
}

head(clusters[[1]])

#6 Combine clusters and tissue info 
bcs <- list()

for (i in 1:length(sample_names)) {
  bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  bcs[[i]]$imagerow <- bcs[[i]]$imagerow   # scale tissue coordinates for lowres image
  bcs[[i]]$imagecol <- bcs[[i]]$imagecol 
  bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
  bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE)
  bcs[[i]]$height <- height$height[i]
  bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names

head(bcs[[1]])

#7 Read in the matrix, barcodes, and genes (it takes 10 mins)
matrix <- list()

for (i in 1:length(sample_names)) {
  matrix[[i]] <- as.data.frame(t(Read10X_h5(matrix_paths[i])))
}

head(matrix[[1]])

#8 Make summary data.frames

#Total UMI per spot
umi_sum <- list() 

for (i in 1:length(sample_names)) {
  umi_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                             sum_umi = Matrix::rowSums(matrix[[i]]))
  
}
names(umi_sum) <- sample_names

umi_sum <- bind_rows(umi_sum, .id = "sample")
head(umi_sum)

#Total Genes per spot
gene_sum <- list() 

for (i in 1:length(sample_names)) {
  gene_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                              sum_gene = Matrix::rowSums(matrix[[i]] != 0))
  
}
names(gene_sum) <- sample_names

gene_sum <- bind_rows(gene_sum, .id = "sample")
head(gene_sum)

#9 Merge all the necessary data
bcs_merge <- bind_rows(bcs, .id = "sample")
bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))
bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))
head(bcs_merge)

