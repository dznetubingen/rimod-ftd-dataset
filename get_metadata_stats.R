##########
# Collect statistics about mapping and read counts for the sequencing arrays
##########
library(stringr)

setwd("~/dzne/rimod")

###
# CAGE-seq statistics
###
md <- read.table("figures_ftd_dataset/rimod_ftd_dataset_table_filtered.txt", sep="\t", header=T)
mapping <- read.table("RiMod_ID_mapping.txt", sep="\t", header=T)
stats = read.csv("data/FTD_Brain_corrected.csv", stringsAsFactors = F)

md <- merge(md, mapping, by.x="RimodID", by.y="new_id")

# format stats table
stats <- stats[stats$REGION == "frontal",]
stats$SAMPLEID <- str_pad(stats$SAMPLEID, 5, side="left", pad="0")
stats$SAMPLEID[stats$SAMPLEID == "A144_12"] <- "0A144"

md <- merge(md, stats, by.x="old_id", by.y="SAMPLEID")
print("CAGE statistics")
print(mean(as.numeric(na.omit(md$Number_of_input_reads))))
print(mean(as.numeric(na.omit(md$Number_of_Uniquely_mapped_reads))))
print(mean(as.numeric(na.omit(md$Ratio_of_uniquely_mapped_reads))))
print(min(as.numeric(na.omit(md$Ratio_of_uniquely_mapped_reads))))


rin = as.numeric(na.omit(md$RIN))

###
# RNA-seq Quality plots
###

rna <- read.table("data/rnaseq/multiqc_data/multiqc_general_stats.txt", sep="\t", header=T, stringsAsFactors = F)
rna <- rna[!grepl("_2", rna$Sample),]
umapped <- as.numeric(na.omit(rna$STAR_mqc.generalstats.star.uniquely_mapped))
print("RNA-seq")
print(mean(umapped))
print(sd(umapped))
print(min(umapped))
print(mean(rna$STAR_mqc.generalstats.star.uniquely_mapped_percent))
print(min(rna$STAR_mqc.generalstats.star.uniquely_mapped_percent))


###
# smRNA-seq quality statistics
srna1 <- read.table("data/smrnaseq/multiqc_general_stats_2018.txt", sep="\t", header=T, stringsAsFactors = F)
srna2 <- read.table("data/smrnaseq/multiqc_general_stats_2019.txt", sep="\t", header=T, stringsAsFactors = F)
srna = data.frame(sample = c(srna1$Sample, srna2$Sample), 
                  "uniquely_mapped" = c(srna1$STAR_mqc.generalstats.uniquely_mapped, srna2$STAR_mqc.generalstats.star.uniquely_mapped),
                  "total_sequences" = c(srna1$FastQC_mqc.generalstats.total_sequences, srna2$FastQC_mqc.generalstats.fastqc.total_sequences))
total_sequences = as.numeric(na.omit(srna$total_sequences))
# clean up OASIS2 extra files
srna <- srna[!grepl("--", srna$sample),]
srna <- srna[!grepl("-trimmed", srna$sample),]
srna <- srna[!grepl("allspeciesCounts", srna$sample),]
srna <- srna[!grepl("std.err", srna$sample),]
srna <- srna[!grepl("std", srna$sample),]
srna_mean = mean(srna$uniquely_mapped)
srna_sd = sd(srna$uniquely_mapped)
print("smRNAseq")
print(srna_mean)
print(srna_sd)




