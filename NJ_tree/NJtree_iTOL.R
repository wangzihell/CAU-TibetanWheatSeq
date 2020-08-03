#!/usr/bin/env Rscript

# generate file for iTOL

# libraries
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-d", "--distF"), dest = "distF", default = "",
              help = "dist matrix file, n rows by n cols."),
  make_option(c("-s", "--sampleF"), dest = "sampleF", default = "",
              help = "sample name that matches dist file."),
  make_option(c("-a", "--annoF"), dest = "annoF", default = "",
              help = "annotation file, which contains sample name to show, sample group."),
  make_option(c("-p", "--sampleP"), dest = "sampleP", default = "",
              help = "only sample name in this file will be plot."),
  make_option(c("-g", "--tmpG"), dest = "tmpG", default = NULL,
              help = "tmp group, if it exist, groups in this file will replace that in metadata."),
  make_option(c("-c", "--colorF"), dest = "colorF", default = "colorFile_clade.txt",
              help = "group name and color."),
  make_option(c("-O", "--outgroup"), dest = "outgroup", default = "",
              help = "group name and color."),
  make_option(c("-t", "--title"), dest = "titleP", default = NULL,
              help = "group name and color."),
  make_option(c("-P", "--plot"), dest = "plotT", default = F,
              help = "plot tree on pdf."),
  make_option(c("-o", "--out"), dest = "out", default = "output",
              help = "output prefix.")
)

parser <- OptionParser(usage = "Rscript treeWholeNJ.R [options] Tip group files",
                       description = 
                       "
if no dist and sample name files provided, use whole genome data. 
Tip group files is used to plot points on branch ends, if none, no points plotted.
                       ",
                       option_list = option_list)

arguments <- parse_args(parser, positional_arguments=c(0,Inf))

distF <- arguments$options$distF
sampleF <- arguments$options$sampleF
annoF <- arguments$options$annoF
sampleP <- arguments$options$sampleP
colorF <- arguments$options$colorF
tmpG <- arguments$options$tmpG
plotT <- arguments$options$plotT
titleP <- arguments$options$titleP
out <- arguments$options$out

# read in data
dist_mat <- read.table(distF)

SM <- read.table(sampleF, header = F, stringsAsFactors = F)
colnames(dist_mat) = SM[[1]]
rownames(dist_mat) = SM[[1]]

meta_data <- read.table(annoF, header = T, stringsAsFactors = F, sep="\t")
meta_data <- meta_data[complete.cases(meta_data),]

present_sample <- read.table(sampleP, header = F, stringsAsFactors = F)

# pre-process
if(!is.null(tmpG)) {
  tmp_group <- read.table(tmpG, header = F, stringsAsFactors = F)
  meta_data[match(tmp_group[,1], meta_data$Sample),4] <- tmp_group[,2]
}

tmp_index <- rownames(dist_mat) %in% meta_data$Name[match(present_sample[[1]], meta_data$Sample)]
dist_mat <- dist_mat[tmp_index, tmp_index]

sample_present <- meta_data[meta_data$Name %in% rownames(dist_mat),]

LIST <- sample_present$label
names(LIST) <- sample_present$Name
rownames(dist_mat) = LIST[rownames(dist_mat)]

dis_matrix <- as.dist(dist_mat)

# define colors
colorFile <- read.table(colorF, header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
colors <- colorFile[,2]
names(colors) <- colorFile[,1]

# NJ-tree

NJ_tree <- njs(dis_matrix)
outgroup <- arguments$options$outgroup
if (outgroup != ""){
  NJ_tree <- root(NJ_tree, outgroup)
}

sample_present_tree <- sample_present[,c(3,1,2,4)]
DF <- full_join(as_tibble(NJ_tree), sample_present_tree, by = 'label')

if(!plotT){
  tree <- as.treedata(DF)
  write.beast(tree, paste0(out, ".nexus"))

  branchcol <- DF %>% filter(!is.na(label)) %>% mutate(recode(Group, !!!colors))  %>% mutate(type = "branch", line = "normal") %>% select(c(4,9,8,10))

  write("TREE_COLORS\nSEPARATOR TAB\nDATA", paste0(out, ".branchcol.txt"))
  write.table(branchcol, paste0(out, ".branchcol.txt"), append = T, col.names = F, row.names = F, quote = F, sep = "\t")
} else {
  col.vector <- vector(mode="character",length=nrow(NJ_tree$edge))
  n.tips <- length(NJ_tree$tip.label)
  col.vector[NJ_tree$edge[,2]>n.tips] <- "black"
  edge.data <- as.data.frame(NJ_tree$edge)
  for(i in seq_along(NJ_tree$tip.label)){
    edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
    col.vector[edge.row] <- colors[DF$Group[match(NJ_tree$tip.label[i], DF$label)]]
  }
  tip.color <- colors[DF$Group]
  names(tip.color) <- as.character(DF$label)
  pdf( paste0(out, ".pdf"), width = 10, height = 10)
  plot(NJ_tree, type = "u", show.tip.label = T, edge.color = col.vector, cex=0.2)
  legend("bottomleft", fill = colors, legend = names(colors), border = F, box.col = "white", horiz=T, cex = 1.3, x.intersp = 0.3)
  if(!is.null(titleP)) title(titleP, cex.main = 2)
  dev.off()
}
