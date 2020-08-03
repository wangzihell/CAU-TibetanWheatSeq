#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-A", "--groupA"), dest = "groupfA", default = "",
              help = "group A"),
  make_option(c("-B", "--groupB"), dest = "groupfB", default = "",
              help = "group B"),
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "chrom"), 
  make_option(c("-b", "--binsize"), dest = "binsize", default = "100",
              help = "binsize of data, in kb, [default: 100]"), 
  make_option(c("-t", "--title"), dest = "title", default = "",
              help = "[opt] title of plot [default: no title]"), 
  make_option(c("-w", "--wd"), dest = "wd", default = "/data2/rawdata2/readDepth",
              help = "[opt] DP file dir [default: /data2/rawdata2/readDepth]"),
  make_option(c("-s", "--suffix"), dest = "suffix", default = ".100k.norm",
              help = "[opt] DP file suffix. [default: .100k.norm]"),
  make_option(c("-d", "--dataOnly"), dest = "dataOnly", action="store_true", default = FALSE, 
              help="Print little output")
)

parser <- OptionParser(usage = "cgmaptools heatmap [options]",
                       option_list=option_list, description = "      (aka mCBinHeatmap)\
                       Description: Plot methylation dynamics of target region for multiple samples [heatmap]\
                       Contact:     Zhu, Ping; pingzhu.work@gmail.com\
                       Last update: 2017-09-16\
                       Example: \
                       mCBinHeatmap.R -i input -m white -o chr1.xxx-xxx.pdf \
                       -Input File Format: \
                       1st line is the header.\
                       Each column contains methylation measurements of a sample. \
                       Example: \
                       Sample1  Sample2 ...  \
                       0.1      0.1     ...  \
                       0.1      0.1     ...  \
                       "
)

arguments <- parse_args(parser)
CHR <- arguments$chr
WD <- arguments$wd
suffix <- arguments$suffix
binsize <- as.numeric(arguments$binsize)

groupfA <- arguments$groupfA
if(groupfA == ""){ # default, STDIN
  print_help(parser)
} else { # user specified
  if( file.access(groupfA) == -1){ # file not exists
    print_help(parser)
  }
}

groupfB <- arguments$groupfB
if(groupfB == ""){ # default, STDIN
  print_help(parser)
} else { # user specified
  if( file.access(groupfB) == -1){ # file not exists
    print_help(parser)
  }
}

groupA <- read.delim(groupfA, header = F, sep = "\t", stringsAsFactors = F)
groupB <- read.delim(groupfB, header = F, sep = "\t", stringsAsFactors = F)

#
lA <- lapply(gsub("ZY-", "S", unlist(groupA[,1])), function(i) {r <- read.table(paste0(WD, "/", i,"/", CHR, suffix) ,sep="\t", stringsAsFactors = F, header = T)})
groupA <- do.call(cbind , lA)
lB <- lapply(gsub("ZY-", "S", unlist(groupB[,1])), function(i) {r <- read.table(paste0(WD, "/", i,"/", CHR, suffix) ,sep="\t", stringsAsFactors = F, header = T)})
groupB <- do.call(cbind , lB)

CNVindex <- rep(NA, nrow(groupA))
Amean <- rep(NA, nrow(groupA))
Bmean <- rep(NA, nrow(groupA))

# some sample with high DP will influence results
groupA[groupA > 1] <- 1
groupB[groupB > 1] <- 1

for(i in 2:nrow(groupA)) {
  #pval <- abs(log(t.test(groupA[i,],groupB[i,])$p.value, base = 10))
  Amean[i] <- rowMeans(groupA[i,])
  Bmean[i] <- rowMeans(groupB[i,])
  CNVindex[i] <- Amean[i] - Bmean[i]
}
CNVindex <- unlist(CNVindex)

if (!isTRUE(arguments$dataOnly)){
  p <- ggplot() + 
    geom_point(aes(x = 1:length(CNVindex), y = CNVindex), shape = 16, alpha = 0.7, size = 0.7) +
    labs(title = arguments$title, x = "POS") +
    theme_cowplot() +
    scale_x_continuous(breaks = seq(0, 9000, 500), labels = function(x) paste0(x/10, "M"))

  ggsave(paste0(CHR, ".CNV-index.png"), p, height = 9, width = 16)
} else {
  df <- data.frame(CHR = CHR, BIN_START = seq(1, (length(CNVindex)-1) * binsize * 1000 + 1, binsize * 1000), BIN_END = seq(binsize * 1000, length(CNVindex) * binsize * 1000, binsize * 1000), CNVindex = CNVindex, Amean = Amean, Bmean = Bmean)
  
  # some location have same DP
  df[is.na(df)] <- 0
  colnames(df) <- c("CHR", "BIN_START", "BIN_END", "CNVindex", "Amean", "Bmean")
  # turn numeric to integer to turn off scitific notaion in write.csv
  class(df$BIN_START) <- "integer"
  class(df$BIN_END) <- "integer"
  write.table(df, file = paste0(CHR, ".CNV-index.xls"), append = F, quote = F, sep = "\t", row.names = F)
}
