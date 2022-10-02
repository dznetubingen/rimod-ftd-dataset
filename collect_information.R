library(stringr)

#####################
# Collect mutational information
#####################
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

######################
# Format count tables for upload
#######################
mapping <- read.table("~/dzne/rimod/RiMod_ID_mapping.txt", header=T)
save_path = "/home/kevin/dzne/rimod/data/count_tables/"

# format RNA-seq
rna <- read.table("~/dzne/rimod/results/rnaseq/frontal_lengthScaledTPM_counts.txt", sep="\t", header=T, row.names = 1)

cols <- colnames(rna)
cols <- gsub("X", "", cols)
cols <- str_split(cols, pattern="_", simplify = T)[,1]
cols <- str_pad(cols, width=5, side="left", pad="0")
keep <- !grepl("5108", cols)
cols <- cols[keep]
rna <- rna[,keep]
rna.map <- mapping[mapping$old_id %in% cols,]
rna.map <- rna.map[order(match(rna.map$old_id, cols)),]
cols.remapped <- as.character(rna.map$new_id)
colnames(rna) <- cols.remapped

write.table(rna, paste0(save_path, "RiMod_RNAseq_lengthScaledTPM_counts.txt"), sep="\t", quote=F)

# format smRNA-seq sample names
mirna <- read.table("~/dzne/rimod/data/smrnaseq/rimod_human_frontal_smRNAseq_counts.txt", sep="\t", header=T, row.names=1)

cols <- colnames(mirna)
cols <- gsub("X", "", cols)
cols <- gsub("sample_", "", cols)
cols <- str_sub(cols, 1, 5)
mirna.map <- mapping[mapping$old_id %in% cols,]
mirna.map <- mirna.map[order(match(mirna.map$old_id, cols)),]
cols.remapped <- as.character(mirna.map$new_id)
colnames(mirna) <- cols.remapped

write.table(mirna, paste0(save_path, "RiMod_smRNAseq_counts.txt"), sep="\t", quote=F)

# format CAGE-seq
cage <- read.table("~/dzne/rimod/data/cage/Raw_All7RegionSamps_3kbGR_DF.txt", sep="\t", header=T)

rownames(cage) <- cage$Clus_ID
cage <- cage[, grepl("_fro_", colnames(cage))]
cols <- colnames(cage)
cols <- str_split(cols, pattern="_", simplify = T)[,2]
cols <- str_pad(cols, width=5, side="left", pad="0")
cage.map <- mapping[mapping$old_id %in% cols,]
cage.map <- cage.map[order(match(cage.map$old_id, cols)),]
cols.remapped <- as.character(cage.map$new_id)
colnames(cage) <- cols.remapped

write.table(cage, paste0(save_path, "RiMod_CAGEseq_cluster_counts.txt"), sep="\t", quote=F)

# format methylation data
met <- read.table("~/dzne/rimod/data/methylation/rimod_frontal_methylation/mVals_matrix_frontal_methylation.txt", sep="\t", header=T)
cols <- colnames(met)
cols <- gsub("X", "", cols)
cols <- str_split(cols, pattern="_", simplify = T)[,1]
cols <- str_pad(cols, width=5, side="left", pad="0")
cols[cols == "14412"] <- "0A144"
met.map <- mapping[mapping$old_id %in% cols,]
met.map <- met.map[order(match(met.map$old_id, cols)),]
cols.remapped <- as.character(met.map$new_id)
colnames(met) <- cols.remapped

write.table(met, paste0(save_path, "RiMod_methylation_mValues.txt"), sep="\t", quote=F)

