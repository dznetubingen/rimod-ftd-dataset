##################################################
### DESeq2 analysis of CAGE-seq clusters to be used
### with TF-activity pipeline                        
##################################################


##########################
## PARAMETER SECTION

# load libs
library(DESeq2)
library(stringr)
library(viridis)
library(pheatmap)
library(ggplot2)
library(factoextra)
### Hard-coded section
script_name = "cage_cluster_deseq_analysis.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = 10
metadata = "~/dzne/rimod/data/FTD_Brain_corrected.csv"

count_file = "~/dzne/rimod/data/cage/tf_enrichment_analysis_050420/RiMod_CAGEseq_cluster_count_table_frontal.txt"
analysis_dir = "~/dzne/rimod/data/cage/tf_enrichment_analysis_050420/"
region <- "fro"

###########################

# set working directory
setwd(analysis_dir)
# Create sub-folder for current analysis
dir.create(paste("CAGE_deseq_analysis","_",region, "_", current_time, sep=""))
setwd(paste("CAGE_deseq_analysis", "_", region, "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, count_file, analysis_dir, script_name, region)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "Count_file", "Analysis_directory", "Script_name", "Region")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")


# Read cage genewise count table (created by Tenzin)
cage <- read.table(count_file, sep="\t", header=T, row.names = 1, stringsAsFactors = F)
# TODO: think about how to remove this hard-coded part ...
# Remove sample sample_09218_froR
# Keep only desired region
cage <- cage[,grepl(region, colnames(cage))]
# Remove 'froR' samples
cage <- cage[,!grepl("froR", colnames(cage))]


# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# bring counts and md in similar format
cage.samples <- as.character(gsub("sample_","",colnames(cage)))
cage.samples <- as.character(sapply(cage.samples, function(x){strsplit(x, split=paste("_", region, sep=""))[[1]][[1]]}))
md <- md[md$SAMPLEID %in% cage.samples,]
md <- md[match(cage.samples, md$SAMPLEID),]


#### REMOVE ALL SPORADIC CASES ####
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
cage <- cage[,keep]
md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe

# PH
ph <- as.numeric(md$PH)
ph.mean <- mean(na.omit(ph))
ph[is.na(ph)] <- ph.mean
md$PH <- ph

rownames(md) <- colnames(cage)
#===========================================#
# DESeq2 analysis
# Generate DDS object
dds <- DESeqDataSetFromMatrix(cage,
                              colData = md,
                              design = ~ PH +  GENDER + DISEASE.CODE)

# Specify control group
dds$DISEASE.CODE <- relevel(dds$DISEASE.CODE, ref = "control")

# apply prefiltering
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
pval_cut <- 0.05
res.mapt <- results(dds, c("DISEASE.CODE", "FTD_MAPT", "control"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <=pval_cut,]
print(dim(deg.mapt))
### GRN - control
res.grn <- results(dds, c("DISEASE.CODE", "FTD_GRN", "control"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("DISEASE.CODE", "FTD_C9", "control"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
print(dim(deg.c9))

###########
## Save results

write.table(res.mapt, paste("deseq_result_mapt.ndc", "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc",  "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only significant genes for online tools
write.table(rownames(deg.mapt), paste("DEGs_mapt.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(deg.grn), paste("DEGs_grn.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(rownames(deg.c9), paste("DEGs_c9.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, row.names=F)

# Divde in up and down regulated genes
# MAPT
mapt.up <- deg.mapt[deg.mapt$log2FoldChange > 0,]
mapt.down <- deg.mapt[deg.mapt$log2FoldChange < 0,]
write.table(rownames(mapt.up), "DEGs_UP_mapt.ndc.txt", quote=F, row.names=F)
write.table(rownames(mapt.down), "DEGs_Down_mapt.ndc.txt", quote=F, row.names=F)
# GRN
grn.up <- deg.grn[deg.grn$log2FoldChange > 0,]
grn.down <- deg.grn[deg.grn$log2FoldChange < 0,]
write.table(rownames(grn.up), "DEGs_UP_grn.ndc.txt", quote=F, row.names=F)
write.table(rownames(grn.down), "DEGs_Down_grn.ndc.txt", quote=F, row.names=F)
# C9orf72
c9.up <- deg.c9[deg.c9$log2FoldChange > 0,]
c9.down <- deg.c9[deg.c9$log2FoldChange < 0,]
write.table(rownames(c9.up), "DEGs_UP_c9.ndc.txt", quote=F, row.names=F)
write.table(rownames(c9.down), "DEGs_Down_c9.ndc.txt", quote=F, row.names=F)

########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, paste("deseq_normalized_counts", "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

# reg log transformed values
rld <- vst(dds, blind=FALSE)
rld.mat <- assay(rld)
write.table(rld.mat, paste("deseq_vst_values","_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

################################
## Plotting section ############

## PCA
pca <- plotPCA(rld, intgroup = "DISEASE.CODE", ntop=5000, returnData = TRUE)
pca$Group <- gsub("_", "-", pca$group)

pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) 
pca
ggsave("~/dzne/rimod/figures_ftd_dataset/pca/CAGE_PCA_deseq_vct.png", height=6, width=6, dpi=300)

# Calculate the row variance, then sort the matrix by variance and keep the top 5000 genes
# Use the resulting matrix to calculate the PCA
n_top = 5000
mat <- assay(rld)
vars = rowVars(mat)
mat <- mat[sort(vars, index.return = TRUE, decreasing = TRUE)$ix,]
mat <- mat[1:n_top,]

all.pca <- prcomp(t(mat), retx = T)
importance <- as.data.frame(summary(all.pca)$importance)
pc1_pct_variance <- round(importance$PC1[2], digits=2)
pc2_pct_variance <- round(importance$PC2[2], digits=2)

df <- data.frame(PC1 = all.pca$x[,1], PC2 = all.pca$x[,2], Group=as.character(rld$DISEASE.CODE ))
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) + 
  xlab(paste0("PC1 (", pc1_pct_variance, "%)")) +
  ylab(paste0("PC2 (", pc2_pct_variance, "%)"))
pca

ggsave("~/dzne/rimod/figures_ftd_dataset/pca/CAGE_PCA_deseq_vst_v2.png", height=6, width=6, dpi=300)
