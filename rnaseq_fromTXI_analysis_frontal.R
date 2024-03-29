###########################################################
# Analysis of Salmon quantified RiMod frontal RNA-seq data
##########################################################
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(pheatmap)
library(ggplot2)

rimod_base_dir = "/home/kevin/dzne/rimod/"

setwd(paste0(rimod_base_dir, "results/rnaseq/"))

#### Hard-coded section
script_name = "rnaseq_fromXI_analysis_rimod_frontal.R"
date = Sys.Date()
current_time = gsub(":", ".", gsub(" ", "_", Sys.time()))
####

# parameters parsing
row_sum_cutoff = 10
row_sum_samples_nr = 5
metadata = file.path(rimod_base_dir, "data/FTD_Brain_corrected.csv")
analysis_dir = getwd()
region <- "fro"
salmon_files = file.path(rimod_base_dir, "results/rnaseq/frontal_lengthScaledTPM_counts.txt")
#====================================================================#

# Create sub-folder for current analysis
dir.create(paste("RNAseq_analysis","_",region, "_", current_time, sep=""))
setwd(paste("RNAseq_analysis", "_", region, "_",current_time, sep=""))

# Save parameters in config file
params <- c(current_time, as.character(row_sum_cutoff), metadata, salmon_files, analysis_dir, script_name, region)
param.names <- c("Time", "Row_sum_cutoff", "Metadata", "salmon_files", "Analysis_directory", "Script_name", "Region")
params <- data.frame(param.name = param.names, param = params)
write.table(params, paste("config_file", current_time, sep="_"), quote = F, row.names = F, sep="\t")

#=================================================================#

# Load metadata
md <- read.csv(metadata, stringsAsFactors = FALSE)
md$SAMPLEID <- as.character(sapply(md$SAMPLEID, function(x){strsplit(x, split="_")[[1]][[1]]}))
md$SAMPLEID <- str_pad(md$SAMPLEID, width = 5, side = "left", pad = "0") # fill sample ids to 5 digits

# load counts
cts <- read.table(salmon_files, sep="\t", header=T, row.names=1)

# bring counts and md in similar format
rna.samples <- as.character(sapply(colnames(cts), function(x){strsplit(x, split="_")[[1]][[1]]}))
rna.samples <- str_pad(gsub("X", "", rna.samples), width=5, side='left', pad='0')
md <- md[md$SAMPLEID %in% rna.samples,]
md <- md[match(rna.samples, md$SAMPLEID),]

#### REMOVE ALL SPORADIC CASES ####
disease.codes <- c("FTD-C9", "FTD-MAPT", "FTD-GRN", "control")
keep <- md$DISEASE.CODE %in% disease.codes
md <- md[keep,]
# subset TXI
cts <- cts[,keep]

# remove sample 05180 (new Analysis 14.01.2020) --> this sample is a strong outlier in the PCA
keep <- !grepl("5108", colnames(cts))
md <- md[keep,]
cts <- cts[,keep]

md$DISEASE.CODE <- gsub("-", "_", md$DISEASE.CODE) # make disease code names safe
# Split Age covariate into bins
age_bins = 4
md$AGE.BIN <- make.names(cut(md$AGE, breaks=age_bins))

# pmd
md$PMD.MIN. <- as.numeric(md$PMD.MIN.)
pmd.mean <- mean(na.omit(md$PMD.MIN.))
md$PMD.MIN.[is.na(md$PMD.MIN.)] <- pmd.mean
md$pmd <- md$PMD.MIN.

# PH
ph <- as.numeric(md$PH)
ph.mean <- mean(na.omit(ph))
ph[is.na(ph)] <- ph.mean
md$PH <- ph

#===========================================#
# DESeq2 analysis
# Generate DDS object
cts <- round(cts) # round to integer counts
dds <- DESeqDataSetFromMatrix(cts,
                              colData = md,
                              design = ~ PH + GENDER + DISEASE.CODE)


# Specify control group
dds$DISEASE.CODE <- relevel(dds$DISEASE.CODE, ref = "control")

# apply prefiltering
dds <- estimateSizeFactors(dds)
keep <- rowSums((counts(dds, normalized=TRUE) >= row_sum_cutoff)) >= row_sum_samples_nr
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
pval_cut <- 0.05
### MAPT - control
res.mapt <- results(dds, c("DISEASE.CODE", "FTD_MAPT", "control"))
res.mapt <- na.omit(res.mapt)
rownames(res.mapt) <- str_split(rownames(res.mapt), pattern="[.]", simplify = T)[,1]
deg.mapt <- res.mapt[res.mapt$padj <= pval_cut,]
#deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= 0.6,]
print(dim(deg.mapt))

### GRN - control
res.grn <- results(dds, c("DISEASE.CODE", "FTD_GRN", "control"))
res.grn <- na.omit(res.grn)
rownames(res.grn) <- str_split(rownames(res.grn), pattern="[.]", simplify = T)[,1]
deg.grn <- res.grn[res.grn$padj <= pval_cut,]
#deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= 0.6,]
print(dim(deg.grn))

### C9orf72 - control
res.c9 <- results(dds, c("DISEASE.CODE", "FTD_C9", "control"))
res.c9 <- na.omit(res.c9)
rownames(res.c9) <- str_split(rownames(res.c9), pattern="[.]", simplify = T)[,1]
deg.c9 <- res.c9[res.c9$padj <= pval_cut,]
#deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= 0.6,]
print(dim(deg.c9))

###########
## Save results
# Adjust rownames
write.table(res.mapt, paste("deseq_result_mapt.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.grn, paste("deseq_result_grn.ndc",  "_", region, "_", current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(res.c9, paste("deseq_result_c9.ndc", "_", region, "_",current_time, ".txt", sep=""), sep="\t", quote=F, col.names = NA)

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
## Generate count table and vst table
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
pca <- plotPCA(rld, intgroup = "DISEASE.CODE")
png(paste("pca_group_deseq_rLogvals", "_", current_time, ".png", sep=""), width = 1200, height = 900)
pca
dev.off()

pca <- plotPCA(rld, intgroup = "GENDER")
png(paste("pca_gender_deseq_rLogvals", "_", current_time, ".png", sep=""), width = 1200, height = 900)
pca
dev.off()

##### HEATMAPs ################

## MAPT
mapt.vst <- rld.mat[rownames(deg.mapt),]
pheatmap(mapt.vst, scale="row")

#==============================#


# Make more PCAs
# separate PCAs for the different mutation
library(factoextra)
dc <- md$DISEASE.CODE
mapt.rld <- rld.mat[,dc %in% c('control', 'FTD_MAPT')]
grn.rld <- rld.mat[,dc %in% c('control', 'FTD_GRN')]
c9.rld <- rld.mat[,dc %in% c('control', 'FTD_C9')]


# MAPT - control
mapt.pca <- prcomp(t(mapt.rld), retx=T)
mapt.dc <- dc[dc %in% c('control', 'FTD_MAPT')]
mapt.gene <- md$GENE[dc %in% c('control', 'FTD_MAPT')]
mapt.gene[mapt.dc == 'control'] <- 'control'
fviz_eig(mapt.pca)
mapt.x <- as.data.frame(mapt.pca$x)
mapt.x$Disease_code <- mapt.gene
mpca <- ggplot(mapt.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_mapt_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
mpca
dev.off()


# GRN - control
grn.pca <- prcomp(t(grn.rld), retx=T)
grn.dc <- dc[dc %in% c('control', 'FTD_GRN')]
fviz_eig(grn.pca)
grn.x <- as.data.frame(grn.pca$x)
grn.x$Disease_code <- grn.dc
gpca <- ggplot(grn.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_grn_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
gpca
dev.off()

# C9 - control
c9.pca <- prcomp(t(c9.rld), retx=T)
c9.dc <- dc[dc %in% c('control', 'FTD_C9')]
fviz_eig(c9.pca)
c9.x <- as.data.frame(c9.pca$x)
c9.x$Disease_code <- c9.dc
cpca <- ggplot(c9.x, aes(x=PC1, y=PC2, color=Disease_code)) +
  geom_point(size=3) +
  stat_ellipse()
png(paste("pca_c9orf72_rlog", "_", current_time, ".png", sep=""), width = 1200, height = 900)
cpca
dev.off()


# PCA for all samples
all.pca <- prcomp(t(rld.mat), retx = T)
fviz_eig(all.pca)
all.x <- as.data.frame(all.pca$x)
all.x$Disease_code <- dc

# Plot the PCA
df = all.x
df$Group = df$Disease_code
dfGroup = gsub("_", "-", df[''])
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) 
pca

plot_dir = "~/dzne/rimod/figures_ftd_dataset/"
ggsave(paste0(plot_dir, "RNAseq_PCA_without_outlier.png"), width=6, height=6, dpi=300)


# Calculate the row variance, then sort the matrix by variance and keep the top 5000 genes
# Use the resulting matrix to calculate the PCA
mat <- assay(rld)
vars = rowVars(mat)
mat <- mat[sort(vars, index.return = TRUE, decreasing = TRUE)$ix,]

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

ggsave("~/dzne/rimod/figures_ftd_dataset/pca/RNAseq_PCA_without_outlier_v2.png", height=6, width=6, dpi=300)
