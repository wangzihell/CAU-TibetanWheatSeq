#!/usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
CHR <- args[1]
df1 <- read.table(paste0(CHR,".pca"), header = T, stringsAsFactors = F)
var_exp <- unlist(df1[1,c(2:11), drop = TRUE])
df1 <- df1[-1,]
# df2 <- tidyr::gather(df1, PC, value, -Sample)

# metadata
metaF <- "/data2/rawdata2/sample_metadata/withXN/metadata_all_v4.txt"
metadata <- read.table(metaF, header = T, stringsAsFactors = F, sep = "\t")
df1$Group <- metadata[match(df1$Sample, metadata$Sample), "Group"]
df1$label <- metadata[match(df1$Sample, metadata$Sample), "label"]

ggplot(df1, aes(PC1, PC2, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_point() +
  labs(title = paste0(CHR, " PC1+PC2")) + xlab(label=paste0("PC1 (", var_exp[1],"%)")) + ylab(label=paste0("PC2 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC12_without_label.pdf"), width = 7, height = 6.5)

ggplot(df1, aes(PC1, PC2, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_text(aes(label = label), col = "gray", size=1, nudge_y=-0.002) +
  geom_point() +
  labs(title = paste0(CHR, " PC1+PC2")) + xlab(label=paste0("PC1 (", var_exp[1],"%)")) + ylab(label=paste0("PC2 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC12_with_label.pdf"), width = 7, height = 6.5)

ggplot(df1, aes(PC2, PC3, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_point() +
  labs(title = paste0(CHR, " PC2+PC3")) + xlab(label=paste0("PC2 (", var_exp[1],"%)")) + ylab(label=paste0("PC3 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC23_without_label.pdf"), width = 7, height = 6.5)

ggplot(df1, aes(PC2, PC3, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_text(aes(label = label), col = "gray", size=1, nudge_y=-0.002) +
  geom_point() +
  labs(title = paste0(CHR, " PC2+PC3")) + xlab(label=paste0("PC2 (", var_exp[1],"%)")) + ylab(label=paste0("PC3 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC23_with_label.pdf"), width = 7, height = 6.5)

ggplot(df1, aes(PC1, PC3, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_point() +
  labs(title = paste0(CHR, " PC1+PC3")) + xlab(label=paste0("PC1 (", var_exp[1],"%)")) + ylab(label=paste0("PC3 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC13_without_label.pdf"), width = 7, height = 6.5)

ggplot(df1, aes(PC1, PC3, color=Group)) +
  scale_color_manual(values = c(NCC = "#33004D", NCL = "#1E89BB", TS = "#168941", TL = "#D51417", CC = "#FF83FA", CL = "#E6B323")) +
  geom_text(aes(label = label), col = "gray", size=1, nudge_y=-0.002) +
  geom_point() +
  labs(title = paste0(CHR, " PC1+PC3")) + xlab(label=paste0("PC1 (", var_exp[1],"%)")) + ylab(label=paste0("PC3 (", var_exp[2],"%)")) +
  theme_bw()
ggsave(paste0(CHR, "_PC13_with_label.pdf"), width = 7, height = 6.5)
