# Collect mutational information
library(stringr)
mapping <- read.table("~/dzne/rimod/RiMod_ID_mapping.txt", header=T)
md <- read.table("~/dzne/rimod/rimod_ftd_dataset_table.txt", header=T)
mut <- read.csv("~/dzne/rimod/data/FTD_Brain_corrected.csv")

md <- merge(md, mapping, by.x="RimodID", by.y="new_id")

samples <- as.character(mut$SAMPLEID)
samples <- str_pad(samples, width=5, side="left", pad="0")
samples[samples == "A144_12"] <- "0A144"
mut$id <- samples
mut <- mut[, colnames(mut) %in% c("id", "GENE")]

md <- merge(md, mut, by.x="old_id", by.y="id", all.x=TRUE)
md <- md[!duplicated(md$RimodID),]
md <- md[, -1]
colnames(md)[ncol(md)] <- "Mutation"

write.table(md, "~/dzne/rimod/rimod_ftd_dataset_table_v2.txt", sep="\t", row.names = F, quote=F)
