# Generate plots for metadata variables

library(ggplot2)

# load data and drop sporadic samples
md <- read.table("~/dzne/rimod/rimod_ftd_dataset_table.txt", header=T)
groups <- c("control", "FTD-C9", "FTD-GRN", "FTD-MAPT")
md <- md[md$DiseaseCode %in% groups,]
md$Group <- md$DiseaseCode

out_path = "/home/kevin/dzne/rimod/figures_ftd_dataset/qc_plots"

# Age boxplot
p <- ggplot(md, aes(x=Group, y=Age, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  xlab("Group") +
  theme(legend.position="none")
p
ggsave(paste0(out_path, "/age_boxplot.png"), height=4, width=4, dpi=300)

# PH boxplot
md$PH <- as.numeric(as.character(md$PH))
p <- ggplot(md, aes(x=Group, y=PH, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  xlab("Group") +
  ylab("pH") +
  theme(legend.position="none")
p
ggsave(paste0(out_path, "/ph_boxplot.png"), height=4, width=4, dpi=300)

# RIN boxplot
md$RIN <- as.numeric(as.character(md$RIN))
p <- ggplot(md, aes(x=Group, y=RIN, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  xlab("Group") +
  theme(legend.position="none")
p
ggsave(paste0(out_path, "/rin_boxplot.png"), height=4, width=4, dpi=300)

# PMD boxplot
md$PMD <- as.numeric(as.character(md$PMD))
p <- ggplot(md, aes(x=Group, y=PMD, color=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  xlab("Group") +
  theme(legend.position="none")
p
ggsave(paste0(out_path, "/pmd_boxplot.png"), height=4, width=4, dpi=300)

# Gender boxplot
df <- data.frame(Group="test", Gender="F", Value = 1)
for (g in groups) {
  tmp <- md[md$Group == g,]
  gender <- data.frame(table(tmp$Gender))
  colnames(gender) <- c("Gender", "Value")
  gender$Group <- g
  df <- rbind(df, gender)
}
df <- df[-1,]
p <- ggplot(df, aes(x=Group, y=Value, fill=Gender)) + 
  geom_bar(stat="identity", position="fill") +
  theme_minimal() + 
  ylab("Percentage")
p
ggsave(paste0(out_path, "/gender_barplot.png"), height=4, width=4, dpi=300)

