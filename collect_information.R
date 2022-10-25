library(stringr)

# Utility functions
split_and_pad <- function(x){
  x <- str_split(x, pattern="_", simplify = T)[,1]
  x <- str_pad(x, width=5, side="left", pad="0")
  return(x)
}

#####################
# Collect mutational information
#####################
mapping <- read.table("~/dzne/rimod/RiMod_ID_mapping.txt", header=T)
md <- read.table("~/dzne/rimod/rimod_ftd_dataset_table.txt", header=T)
mut <- read.csv("~/dzne/rimod/data/FTD_Brain_corrected.csv")

md <- merge(md, mapping, by.x="RimodID", by.y="new_id")


# Remove spoardic samples
md <- md[md$DiseaseCode %in% c("control", "FTD-C9", "FTD-MAPT", "FTD-GRN"),]
rimod_samples <- as.character(md$RimodID)

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

# filter sporadic samples
rna <- rna[, colnames(rna) %in% rimod_samples]

write.table(rna, paste0(save_path, "RiMod_RNAseq_lengthScaledTPM_counts.txt"), sep="\t", quote=F)

# format smRNA-seq sample names
mirna <- read.table("~/dzne/rimod/data/smrnaseq/rimod_human_frontal_smRNAseq_counts.txt", sep="\t", header=T, row.names=1)
test <- read.table("~/dzne/rimod/data/smrnaseq/rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", header=T)
cols <- colnames(mirna)
cols <- gsub("X", "", cols)
cols <- gsub("sample_", "", cols)
cols <- str_sub(cols, 1, 5)
mirna.map <- mapping[mapping$old_id %in% cols,]
mirna.map <- mirna.map[order(match(mirna.map$old_id, cols)),]
cols.remapped <- as.character(mirna.map$new_id)
colnames(mirna) <- cols.remapped

# remove sporadic
mirna <- mirna[, colnames(mirna) %in% rimod_samples]

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

# remove sporadic samples
cage <- cage[, colnames(cage) %in% rimod_samples]

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

# remove spoardic samples
met <- met[, colnames(met) %in% rimod_samples]

write.table(met, paste0(save_path, "RiMod_methylation_mValues.txt"), sep="\t", quote=F)


###############
# Collect mapping QC information
###############

# RNA-seq QC stats
rna.qc <- read.table("~/dzne/rimod/data/rnaseq/multiqc_data/multiqc_general_stats.txt", sep="\t", header=T)
rna.salmon <- read.table("~/dzne/rimod/data/rnaseq/results_salmon/multiqc_data/multiqc_general_stats.txt", sep="\t", header=T)
rna.star <- read.table("~/dzne/rimod/data/rnaseq/multiqc_data/multiqc_star.txt", sep="\t", header=T)

# remove duplicates
rna.qc$Sample <- split_and_pad(as.character(rna.qc$Sample))
rna.qc <- rna.qc[!duplicated(rna.qc$Sample),]

rna.star$Sample <- gsub("_1", "", rna.star$Sample)
rna.star$Sample <- split_and_pad(as.character(rna.star$Sample))
rna.star <- rna.star[, c(1,2)]
colnames(rna.star)[2] <- "RNAseq_total_sequences"
rna.star$Sample[rna.star$Sample == "A1442"] <- "0A144"

rna.salmon$Sample <- split_and_pad(rna.salmon$Sample)
rna.stats = merge(rna.qc, rna.salmon, by.x = "Sample", by.y="Sample")
rna.stats <- rna.stats[, c(1, 5, 6, 7, 8, 9)]
colnames(rna.stats) <- c("Sample", 
                         "RNAseq_STAR_uniquely_mapped_percent", 
                         "RNAseq_STAR_uniquely_mapped", 
                         "RNAseq_percent_trimmed", 
                         "RNAseq_Salmon_num_mapped", 
                         "RNAseq_Salomn_percent_mapped")

rna.stats = merge(rna.stats, rna.star, on="Sample")

# smRNA-seq QC stats
mirna1 <- read.table("~/dzne/rimod/data/smrnaseq/multiqc_general_stats_2018.txt", sep="\t", header=T)
mirna2 <- read.table("~/dzne/rimod/data/smrnaseq/multiqc_general_stats_2019.txt", sep="\t", header=T)
colnames(mirna2) <- gsub(".star", "", gsub(".fastqc", "", gsub(".bowtie_1", "", gsub(".htseq_count", "", colnames(mirna2)))))
mirna2 <- mirna2[, order(match(colnames(mirna2), colnames(mirna1)))]
mirna <- rbind(mirna1, mirna2)
# remove unused rows
mirna <- mirna[!grepl("--mirna", mirna$Sample),]
mirna <- mirna[!grepl("--srna", mirna$Sample),]
mirna <- mirna[!grepl("-trimmed", mirna$Sample),]
mirna <- mirna[!grepl("_allspeciesCounts", mirna$Sample),]
mirna <- mirna[!grepl("std", mirna$Sample),]
# format sample names
samples <- as.character(mirna$Sample)
samples <- gsub("RNAomeTb", "", gsub("final_5bp_trimmed_sample_", "", samples))
samples[samples == "110140NDCGFM_mm_smallrna_sr_Farah_D_3"] <- "11040"
samples <- str_sub(samples, 1, 5)
mirna$Sample <- samples
mirna <- mirna[mirna$Sample %in% mapping$old_id,]
mirna <- mirna[, c(1, 2, 3)]
colnames(mirna) <- c("Sample",
                     "smRNAseq_STAR_uniquely_mapped_percent",
                     "smRNAseq_STAR_uniquely_mapped")

# CAGE-seq QC stats
cage <- read.csv("~/dzne/rimod/data/FTD_Brain_corrected.csv")
cage <- cage[cage$REGION == "frontal",]
cage <- cage[, c(4, 11, 12, 13)]
samples <- as.character(cage$SAMPLEID)
samples[samples == "A144_12"] <- "0A144"
samples <- str_pad(samples, width=5, side="left", pad="0")
cage$SAMPLEID <- samples
cage <- cage[cage$SAMPLEID %in% mapping$old_id,]
colnames(cage) <- c("Sample", "CAGEseq_total_sequences", "CAGEseq_STAR_uniquely_mapped", "CAGEseq_STAR_uniquely_mapped_percent")

# Collection methylation QC info
met <- read.table("~/dzne/rimod/data/methylation/detection_pvalues.txt", sep="\t", header=T)

# Merge all QC tables together
qc <- merge(rna.stats, mirna, by="Sample", all=TRUE)
qc <- merge(qc, cage, on="Sample", all=TRUE)

qc <- merge(qc, mapping, by.x="Sample", by.y="old_id")
rownames(qc) <- qc$new_id
qc$Sample <- qc$new_id
qc <- qc[, c(-13, -14)]

# remove sporadic samples
qc <- qc[rownames(qc) %in% rimod_samples,]

write.table(qc, "~/dzne/rimod/rimod_ftd_sequencing_qc_statistics.txt", sep="\t", row.names = F, quote=F)


