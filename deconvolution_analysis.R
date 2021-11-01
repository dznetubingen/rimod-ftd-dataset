################################
# Analysis of deconvolution results
################################
library(stringr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(viridis)

# switch working directory
rimod_base_dir = "/home/kevin/dzne/rimod/"
setwd(paste0(rimod_base_dir, "results/rnaseq/deconvolution/deconvolution_analysis"))

plot_dir = "~/dzne/rimod/figures_ftd_dataset/"
# Setup color palette
# color palette one extra MAPT group
mypal_maptSplit <- c("#67e08a", "#19943d", "#db6e1a", "7570B3")
# Color palette for all groups incl. NDC
mypal <- c()
# color palette only for disease groups
mypal_short <- c("#67e08a", "#db6e1a", "7570B3")


####
# Cell composition plot
####

fracs <- read.table("../cdn_predictions.txt", sep="\t", header=T)
colnames(fracs)[1] <- "sample"
fracs$sample <- gsub("X", "", fracs$sample)
fracs$sample <- str_split(fracs$sample, pattern="_", simplify = T)[,1]

# save for other analysis
fracs_raw <- fracs # save for later analysis (supplements)
fracs_table <- fracs
rownames(fracs_table) <- paste0("sample_", fracs_table$sample)
fracs_table <- fracs_table[,-1]
write.table(fracs_table, "deconvolution_predictions_formatted.txt", sep="\t", quote=F, col.names = NA)

# get design matrix
md <- read.csv("../../../../data/FTD_Brain_corrected.csv")

md <- md[md$REGION == "frontal",]
md$sample <- str_split(md$GIVENSAMPLENAME, pattern="_", simplify = T)[,1]
md <- md[md$sample %in% fracs$sample,]
md <- md[match(fracs$sample, md$sample),]
fracs$group <- as.character(md$DISEASE.CODE)
#fracs$group[as.character(md$GENE) == "P301L"] <- "MAPT-P301L"

# remove sporadic
fracs <- fracs[!fracs$group == "Sporadic-TDP",]

# remove Unknown
fracs <- fracs[,-1]

####
# Calculate average increase/decrase
####
mean_fun <- median
fracs <- fracs[, -1]
mapt <- fracs[fracs$group == "FTD-MAPT",]
grn <- fracs[fracs$group == "FTD-GRN",]
c9 <- fracs[fracs$group == "FTD-C9",]
#mp3 <- fracs[fracs$group == "MAPT-P301L",]

control <- fracs[fracs$group == "control",]
mapt <- mapt[, -ncol(mapt)]
mapt <- apply(mapt, 2, mean_fun)
#mp3 <- mp3[, -ncol(mp3)]
#mp3 <- apply(mp3, 2, mean_fun)
grn <- grn[, -ncol(grn)]
grn <- apply(grn, 2, mean_fun)
control <- control[, -ncol(control)]
control <- apply(control, 2, mean_fun)
c9 <- c9[, -ncol(c9)]
c9 <- apply(c9, 2, mean_fun)




# MAPT
mapt_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- mapt[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  mapt_list <- c(mapt_list, pct)
}

# MAPT-P301L
# mp3_list <- c()
# for (i in 1:length(control)) {
#   ct.mean <- control[i]
#   diff <- mp3[i] - ct.mean
#   pct <- diff / ct.mean
#   print(pct)
#   mp3_list <- c(mp3_list, pct)
# }


# GRN
grn_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- grn[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  grn_list <- c(grn_list, pct)
}

# C9ORF72
c9_list <- c()
for (i in 1:length(control)) {
  ct.mean <- control[i]
  diff <- c9[i] - ct.mean
  pct <- diff / ct.mean
  print(pct)
  c9_list <- c(c9_list, pct)
}

# make data frame
grn <- t(data.frame(grn_list))
mapt <- t(data.frame(mapt_list))
#mp3 <- t(data.frame(mp3_list))
c9 <- t(data.frame(c9_list))
df <- data.frame(rbind(mapt, grn, c9))
df$group <- c("MAPT", "GRN", "C9ORF72")
df <- melt(df)

# Generate the plots
#mypal_maptSplit <- c("#7570B3", "#db6e1a", "#67e08a", "#19943d")
mypal <- c("#7570B3", "#db6e1a", "#19943d")
ggplot(data=df, aes(x=variable, y=value, fill=group)) +
  geom_bar(stat='identity', position = position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 12)) +
  scale_fill_manual(values = mypal) + 
  xlab("") + 
  ylab("Percentage difference to control") + 
  ggtitle("Cell Composition Change in Disease Groups")
ggsave(paste0(plot_dir, "cell_composition_percentage_change.png"), width=8, height=6)


# Regress neuronal fractions with other cell types
cells <- fracs
cells$Neurons <- cells$ExNeurons + cells$InNeurons
cells <- cells[, c(-1, -7)]
cells <- melt(cells, id.vars = c("group", "Neurons"))
colnames(cells) <- c("Group", "Neurons", "Celltype", "Fraction")

ggplot(cells, aes(x=Neurons, y=Fraction, color=Celltype)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  theme_minimal() +
  facet_wrap(~Celltype, nrow=1)
#================================#

#########
# Calculate significance in changes
#########
mapt <- fracs[fracs$group %in% c("FTD-MAPT","MAPT-P301L"),]
grn <- fracs[fracs$group == "FTD-GRN",]
c9 <- fracs[fracs$group == "FTD-C9",]
control <- fracs[fracs$group == "control",]

mapt_pvals <- c()
grn_pvals <- c()
c9_pvals <- c()

celltypes <- colnames(mapt[, -ncol(mapt)])

for (ct in celltypes) {
  print(ct)
  
  cont <- control[, colnames(control) == ct]
  tmp.mapt <- mapt[, colnames(mapt) == ct]
  tmp.grn <- grn[, colnames(grn) == ct]
  tmp.c9 <- c9[, colnames(c9) == ct]
  
  res.mapt <- t.test(cont, tmp.mapt)$p.value
  res.grn <- t.test(cont, tmp.grn)$p.value
  res.c9 <- t.test(cont, tmp.c9)$p.value
  
  mapt_pvals <- c(mapt_pvals, res.mapt)
  grn_pvals <- c(grn_pvals, res.grn)
  c9_pvals <- c(c9_pvals, res.c9)
  
}

pval.df <- data.frame("FTD-MAPT" = mapt_pvals, "FTD-GRN" = grn_pvals, "FTD-C9orf72" = c9_pvals, "Celltype" = celltypes)
write.table(pval.df, "~/rimod/RNAseq/analysis/deconvolution/cell_type_difference_testing.txt", sep="\t", quote=F, row.names = F)

# Generate heatmap of p-values
ct <- pval.df$Celltype
df_heatmap <- pval.df[,-4]
df_heatmap[df_heatmap > 0.05] <- NA
anno_df <- data.frame(ct = ct)
rownames(df_heatmap) <- ct
pheatmap(df_heatmap, color = viridis(200, option="A"), cluster_rows = F, 
         cluster_cols = F, width=6, height=6, angle="45")


###########
# Plot for the supplements
###########

# load data

md <- md[match(fracs_raw$sample, md$sample),]
fracs_raw$Group = md$DISEASE.CODE
fracs_raw$Sample <- fracs_raw$sample
fracs_raw <- fracs_raw[,-1]

df <- melt(fracs_raw, id.vars = c("Group", "Sample"))

colnames(df) <- c("Group", "Sample", "Celltype", "Fraction")

df$Group <- factor(as.character(df$Group), levels=c("control", "FTD-C9", "FTD-MAPT", "FTD-GRN"))

# remove 10166
df <- df[!df$Sample == "10166",]

ggplot(df, aes(fill=Celltype, y=Fraction, x=Sample)) + 
  geom_bar(stat="identity", position="fill") + 
  theme_minimal() + 
  theme(axis.text.x=element_blank()) +
  facet_grid(~ Group, scales="free", space="free_x")

ggsave(paste0(plot_dir, "cell_composition_stacked_barplot.png"), width=8, height=6)
#===================================================================#