# estraggo le coordinate del sito di iniezione

#> which(colData(spe_c)$inj_site == "control")
#[1]   11   24   73   92  134  140  154  160  167  196  206  248  267  290  291  311  343  373  379  387  427  502
#[23]  607  622  653  675  680  704  733  741  780  824  855  896  924  960  985 1033 1066 1069 1124 1153 1155 1168
#[45] 1182 1198 1219 1245 1290 1293 1302 1332 1392 1404 1432 1439 1477 1485 1493 1495 1499 1508 1574 1586 1599 1620
#[67] 1646 1654 1673 1693 1789 1812 1818 1821 1822 1833 1851 1869 1878 1901 1911 1917 1929 1936 1972 1974 1988 2040
#[89] 2096 2175 2185 2208 2224 2236 2249 2263 2294 2332 2348 2368 2408 2427 2442 2526 2545 2552 2553 2599

library(readr)
tissue_positions_list <- read_csv("control/outs/spatial/tissue_positions_list.csv")

n <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
names(tissue_positions_list) <- n

# remove spot not detected (in tissue == 0)
tissue_positions_list <- subset(tissue_positions_list, in_tissue == 1)

# combine tissue position with injenction site 
control_positions <- merge(tissue_positions_list, control_injection_site, by = "Barcode")

# Extract only injection site positions 
inj_positions <- subset(control_positions, injection_site == "control")

# calcolo la distanza utilizzando le coordinate dello spot centrale del sito di iniezione 
# inj_positions[60,]
#                 Barcode in_tissue arrary_row array_col pixel_row pixel_col injection_site
# 1495 GATCTTGGAGGGCATA-1         1         44        70      4424      3995        control

# x = 44 ; y = 70

tissue_positions_list$distance <- NA
i=1
while(i <= 2679){
  tissue_positions_list[i,7] <- sqrt((tissue_positions_list[i,3] - 44)^2 + (tissue_positions_list[i,4] - 70)^2)
  i = i + 1
  }
  

colData(spe_c)$inj_site_distance <- tissue_positions_list$distance

