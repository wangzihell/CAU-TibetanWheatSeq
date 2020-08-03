#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(pophelper))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-d", "--workD"), dest = "workD", default = "",
              help = "directory containing STRUCTURE output Q files"),
  make_option(c("-s", "--sampleF"), dest = "sampleF", default = "",
              help = "sample names matching STRUCTURE output Q files"),
  make_option(c("-a", "--annoF"), dest = "annoF", default = "",
              help = "annotation file, which contains sample name to show, sample group."),
  make_option(c("-b", "--baseN"), dest = "baseN", default = "",
              help = "basename of STRUCTURE output Q files"),
  make_option(c("-c","--color"), dest = "color", default = "",
              help = "[opt] A vector of colours for clusters, sep by comma."),
  make_option(c("-S","--share"), dest = "sharelab", default = T,
              help = "[opt] share label between Ks."),
  make_option(c("-p", "--prefix"), dest = "prefix", default = "",
              help = "[opt] output file name. [default: Structure+systime]"),
  make_option(c("-t","--title"), dest = "figure.title", default = "",
              help = "[opt] plot title. [default: Structure plot of `sample-name`]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 25,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 4,
              help = "[opt] height of figure (inch). [default: 7]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: pdf. [default: png]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 300,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 300]")
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

## check arguments
arguments <- parse_args(parser)
#
workD <- arguments$workD
if(workD == ""){ # default, STDIN
  print_help(parser)
}
#
baseN <- arguments$baseN
if(baseN == ""){ # default, STDIN
  print_help(parser)
}
#
sampleF <- arguments$sampleF
if(sampleF == ""){ # default, STDIN
  print_help(parser)
} else { # user specified
  if( file.access(sampleF) == -1){ # file not exists
    print_help(parser)
  }
}
#
annoF <- arguments$annoF
if(annoF == ""){ # default, STDIN
  print_help(parser)
} else { # user specified
  if( file.access(annoF) == -1){ # file not exists
    print_help(parser)
  }
}
#
colors <- arguments$color
if(colors == ""){ # default, STDIN
  #colors <- c("#1E89BB", "#E6B323", "#168941", "#D51417","#846F2C", "#33004D")
  colors <-  c("#0099E6", "#846F2C", "#FF6600", "#33004D", "#0D660D", "#991A4D", "#FEB1B1", "#FF004D", "#339933", "#FE9800")
}else { 
  colors <- strsplit(arguments$color, split = ",", fixed = T)[[1]]
}
#
prefix = arguments$prefix
if(prefix == ""){ # default, "FragRegView.Date"
  prefix = paste("Structure", Sys.Date(), sep=".")
} else { # user specified
  prefix <- gsub(".pdf$|.png$", "", prefix, perl=T)
}
#
figure.title = arguments$figure.title
if(figure.title == ""){ # default, STDIN
  figure.title = paste0("STRUCTURE plot of ", baseN)
} 
#
sharelab <- arguments$sharelab
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- arguments$figure.format
figure.resolution <- arguments$figure.resolution
#
if(! figure.format %in% c("pdf", "png")){ # format not support
  print_help(parser)
} else {
  outfile <- prefix
}
#
qfile <- readQ(paste(workD, list.files(path = workD, pattern = baseN, recursive = T), sep = "/"))
sampleN <- read.table(sampleF, stringsAsFactors = F)
metadata <- read.table(annoF, stringsAsFactors = F, sep = "\t", header = T)
#
sample_present <- metadata[match(sampleN[[1]], metadata$Name),]
#
tmp_index <- order(sample_present$Group)
sample_present  <- sample_present[tmp_index,]
#
sorted_popfile <- lapply(qfile, function(df){
  df <- df[tmp_index,]
  rownames(df) <- sample_present$label
  df
})
#
fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(sorted_popfile,fn1))
sorted_popfile <- alignK(sorted_popfile)

#
if(sharelab){
  plotQ(sorted_popfile, 
        imgoutput="join", 
        showindlab=T, showlegend = T, useindlab=T, showtitle = T, showsubtitle = T,
        grplab = data.frame(Group = sample_present[,4], stringsAsFactors = F), 
        splab = spnames, splabsize = 15, titlesize = 13, subtitlesize = 13,
        clustercol = colors, sortind = NA, sharedindlab = T,
        height = figure.height, width = figure.width, imgtype = figure.format, outputfilename = outfile,
        titlelab = figure.title, subtitlelab = "Downsample 10% first; filter criteria: MAF>0.05; Missing-rate<0.15; LD(R^2)<0.5; ",
        grplabsize = 5, legendtextsize = 5, divsize = 0.5, divtype = 1, divcol = "black",
        indlabheight=0.18, indlabspacer=1, indlabsize = 2)
} else {
  plotQ(sorted_popfile, 
        imgoutput="join", 
        showindlab=T, showlegend = T, useindlab=T, showtitle = T, showsubtitle = T,
        grplab = data.frame(Group = sample_present[,4], stringsAsFactors = F), 
        clustercol = colors, sortind = "all", sharedindlab = F,
        height = figure.height, width = figure.width, imgtype = figure.format, outputfilename = outfile,
        titlelab = figure.title, subtitlelab = "thin101 means randomly downsample 10% of SNPs, rep1; thin102, 10%, rep2.",
        grplabsize = 3, legendtextsize = 5, divsize = 1,
        indlabheight=0.18, indlabspacer=1, indlabsize = 2, barsize = 1, barbordersize=0)
}
