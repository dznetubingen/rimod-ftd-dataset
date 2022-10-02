#######################################################
# Analysis of fronal smRNA-seq data of RiMOD project
#######################################################
library(DESeq2)
library(stringr)
library(pheatmap)
library(viridis)
library(limma)
library(ggplot2)

# Parameters
row_sum_cutoff = 5
pval_cutoff = 0.05
lfc_cutoff = 0.6

setwd("~/dzne/rimod/results/smrnaseq/")

plot_dir = "~/dzne/rimod/figures_ftd_dataset/"

# Load Count data
counts <- read.table("~/dzne/rimod/data/smrnaseq/rimod_human_frontal_smRNAseq_counts.txt", sep="\t", header=T, row.names = 1, check.names = F)
# only keep human miRNAs
keep <- grepl("hsa-", rownames(counts))
counts <- counts[keep,]

# Load Metadata
md <- read.table("~/dzne/rimod/data/smrnaseq/rimod_human_frontal_smRNAseq_metadata.txt", sep="\t", header=T, check.names=F, row.names = 1)
md$id <- as.factor(md$id)

# PH
ph <- read.csv("~/dzne/rimod/data/FTD_Brain_corrected.csv", stringsAsFactors = F)
ph <- ph[ph$REGION == "frontal",]
ph <- ph[ph$SAMPLEID %in% md$id,]
ph <- ph[match(md$id, ph$SAMPLEID),]
md$ph <- as.numeric(ph$PH)
ph.mean <- mean(na.omit(md$ph))
md$ph[is.na(md$ph)] <- ph.mean


dds <- DESeqDataSetFromMatrix(counts,
                              colData = md,
                              design = ~ ph + batch +  gender + dc)



# Specify control group
dds$dc <- relevel(dds$dc, ref = "NDC")
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)
resnames <- resultsNames(dds)

#== Extract results ==#
### MAPT - control
res.mapt <- results(dds, c("dc", "FTD.MAPT", "NDC"))
res.mapt <- na.omit(res.mapt)
deg.mapt <- res.mapt[res.mapt$padj <= pval_cutoff,]
deg.mapt <- deg.mapt[abs(deg.mapt$log2FoldChange) >= lfc_cutoff,]

### GRN - control
res.grn <- results(dds, c("dc", "FTD.GRN", "NDC"))
res.grn <- na.omit(res.grn)
deg.grn <- res.grn[res.grn$padj <= pval_cutoff,]
deg.grn <- deg.grn[abs(deg.grn$log2FoldChange) >= lfc_cutoff,]

### C9orf72 - control
res.c9 <- results(dds, c("dc", "FTD.C9", "NDC"))
res.c9 <- na.omit(res.c9)
deg.c9 <- res.c9[res.c9$padj <= pval_cutoff,]
deg.c9 <- deg.c9[abs(deg.c9$log2FoldChange) >= lfc_cutoff,]


# Save all results
write.table(res.mapt, "deseq_result_mapt.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.grn, "deseq_result_grn.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)
write.table(res.c9, "deseq_result_c9.ndc_frontal_smRNAseq.txt", sep="\t", quote=F, col.names = NA)

# Save differentially expressed miRNAs according to specified cutoff
write.table(deg.mapt, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_mapt.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.grn, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_result_grn.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)
write.table(deg.c9, paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"result_c9.ndc_frontal_smRNAseq.txt", sep=""), sep="\t", quote=F, col.names = NA)

# Save only DEGs (without ohter info) for use in Pathway tools
write.table(rownames(deg.mapt), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_mapt.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.grn), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_grn.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)
write.table(rownames(deg.c9), paste("DEGs_P",pval_cutoff,"_LFC",lfc_cutoff,"_c9.ndc_frontal_smRNAseq_miRNAs.txt", sep=""), sep="\t", quote=F, row.names = FALSE)



########################################
## Generate count table and rLog table
########################################

# normalized count values
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "deseq_normalized_counts_temporal_smRNA.txt", sep="\t", quote=F, col.names = NA)

# reg log transformed values
vst.vals <- varianceStabilizingTransformation(dds, blind=FALSE)
vst.mat <- assay(vst.vals)
write.table(vst.mat, "deseq_vst_values_frontal_smRNA.txt", sep="\t", quote=F, col.names = NA)


## PCA
pca <- plotPCA(vst.vals, intgroup = "dc")
png(paste0(plot_dir, "PCA_rimod_frontal_VST_group.png"), width=800, height=600)
pca
dev.off()
pca
plotPCA(rld, intgroup = "dc")

# remove batch effect with limma
design <- model.matrix(~ md$dc)
x_noBatch <- removeBatchEffect(vst.mat, batch = md$batch, design=design)
nb <- vst.vals
assay(nb) <- x_noBatch

plotPCA(nb, intgroup = "dc")
png(paste0(plot_dir, "PCA_sRNA_rimod_frontal_vst_batchCorrected.png"), width=800, height=600)
df <- plotPCA(nb, intgroup = "dc", returnData = TRUE)
dev.off()

df$Group <- gsub("[.]", "-", df$group)
df$Group[df$Group == "NDC"] <- "control"
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) 
pca

ggsave("smRNAseq_PCA_VST_batchCorrected.png", width=6, height=6, dpi=300)

###
# Create new PCA plot with variance explained in the axis labels
mat <- assay(nb)
vars = rowVars(mat)
mat <- mat[sort(vars, index.return = TRUE, decreasing = TRUE)$ix,]

group <- as.character(nb$dc)
group <- gsub("[.]", "-", group)
group[group == "NDC"] <- "control"
group <- factor(group, levels=c("control", "FTD-C9", "FTD-MAPT", "FTD-GRN"))

all.pca <- prcomp(t(mat), retx = T)
importance <- as.data.frame(summary(all.pca)$importance)
pc1_pct_variance <- round(importance$PC1[2], digits=2)
pc2_pct_variance <- round(importance$PC2[2], digits=2)

df <- data.frame(PC1 = all.pca$x[,1], PC2 = all.pca$x[,2], Group=group)
pca <- ggplot(df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) + 
  xlab(paste0("PC1 (", pc1_pct_variance, "%)")) +
  ylab(paste0("PC2 (", pc2_pct_variance, "%)"))
pca

ggsave("~/dzne/rimod/figures_ftd_dataset/pca/smRNAseq_PCA_VST_batchCorrected_v2.png", height=6, width=6, dpi=300)
