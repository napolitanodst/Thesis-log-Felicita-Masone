# --- CALCOLO DELLE DISTANZE per il controllo e heme_1000

# --- Control ---

## Importo le coordinate relative al campione control
library(readr)
tissue_positions_list <- read_csv("control/outs/spatial/tissue_positions_list.csv", 
                                  col_names = FALSE)
## Rinomino le colonne 
n <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
names(tissue_positions_list_control) <- n

## remove spot not detected (in tissue == 0)
tissue_positions_list <- subset(tissue_positions_list_control, in_tissue == 1)

## Importo le informazioni riguardo il sito di iniezione
control_injection_site <- read_csv("control/control_injection_site.csv")

## combine tissue position with injenction site 
control_positions <- merge(tissue_positions_list_control, control_injection_site, by = "Barcode")

## Extract only injection site positions 
inj_positions_control <- subset(control_positions, injection_site == "control")

# calcolo la distanza utilizzando le coordinate dello spot centrale del sito di iniezione 

# inj_positions[60,]
#                 Barcode in_tissue arrary_row array_col pixel_row pixel_col injection_site
# 1495 GATCTTGGAGGGCATA-1         1         44        70      4424      3995        control

# x = 44 ; y = 70

control_positions$distance <- NA
i=1
while(i <= length(control_positions$distance)){
  control_positions[i,"distance"] <- sqrt((control_positions[i,"array_row"] - 44)^2 + (control_positions[i,"array_col"] - 70)^2)
  i = i + 1
}

## Aggiungo le distanze nell'oggetto spe_c, in particolare in colData(spe_c)
colData(spe_c)$inj_site_distance <- control_positions$distance

## Plot 
vis_gene(
  spe = spe_c,
  sampleid = "control",
  geneid = "inj_site_distance", 
  assayname = "logcounts",
  image_id = "lowres",
  alpha = 1,
  point_size = 2
)

# --- heme_1000 ---

## Importo le coordinate relative al campione control
tissue_positions_heme_1000 <- read_csv("heme_1000/outs/spatial/tissue_positions_list.csv", 
                                  col_names = FALSE)
## Rinomino le colonne 
n <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
names(tissue_positions_list_heme_1000) <- n

## remove spot not detected (in tissue == 0)
tissue_positions_list_heme_1000 <- subset(tissue_positions_list_heme_1000, in_tissue == 1)

## Importo le informazioni riguardo il sito di iniezione
heme_1000_injection_site <- read_csv("heme_1000/heme_1000_injection_site.csv")

## combine tissue position with injenction site 
heme_1000_positions <- merge(tissue_positions_list_heme_1000, heme_1000_injection_site, by = "Barcode")

## Extract only injection site positions 
inj_positions_heme_1000 <- subset(heme_1000_positions, injection_site == "heme")

# calcolo le distanze dallo spot [81]
# > inj_positions_heme1000[81,]
#                 Barcode in_tissue arrary_row array_col pixel_row pixel_col injection_site
# 1536 GACGACGATCCGCGTT-1         1         52        70      5118      3867           heme

# x = 52 ; y = 70

heme_1000_positions$distance <- NA
i=1
while(i <= length(heme_1000_positions$distance)){
  heme_1000_positions[i,"distance"] <- sqrt((heme_1000_positions[i,"array_row"] - 52)^2 + (heme_1000_positions[i,"array_col"] - 70)^2)
  i = i + 1
}

## Aggiungo le distanze nell'oggetto spe_c, in particolare in colData(spe_c)
colData(spe_h1000)$inj_site_distance <- heme_1000_positions$distance

## plot
vis_gene(
  spe = spe_h1000,
  sampleid = "heme_1000",
  geneid = "inj_site_distance", 
  assayname = "logcounts",
  image_id = "lowres",
  alpha = 1,
  point_size = 2
)
