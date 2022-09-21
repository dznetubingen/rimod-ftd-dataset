# Collect mutational information

mapping <- read.table("~/dzne/rimod/RiMod_ID_mapping.txt", header=T)
md <- read.table("~/dzne/rimod/rimod_ftd_dataset_table.txt", header=T)
mut <- read.csv("~/dzne/rimod/data/FTD_Brain_corrected.csv")