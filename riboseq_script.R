# NGS / HW5 / Ribo-seq
####################################
# 1. Installing package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

install.packages("tidyverse") 
library("tidyverse") 

# 2. Data retrieval and preprocessing
setwd("/Users/katerinaoleynikova/Downloads")
countdata <- read.csv("counts.csv")
coldata <- read.csv("design.csv")
# remove first column
countdata <- countdata[-1]
# convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# 3. Analysis
# We start the bioconductor DESeq2 analysis worfklow,
# first by initializing the DESeqDataSeet from the matrix we've created above:
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Running the DESeq Pipeline
dds <- DESeq(dds)

# Quality Control Plots
# We next perform some basic QC plots to examine the data behaviour,
# first by looking at dispersions to ensure they behave smoothly and broadly decrease with counts:
options(repr.plot.height=4, repr.plot.width=6)
plotDispEsts(dds, main="Dispersion plot")

# In general we're most interested in logs of the counts; DESeq has a regularlized log transform built in
# Regularized log transformation for clustering/heatmaps, etc:
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Now we'll make sure sure the samples cluster together by condition, first by clustering:
library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(coldata$condition))]
options(repr.plot.height=2, repr.plot.width=6)
DESeq2::plotPCA(rld, intgroup="condition")

# 4. Dif. expression analysis
# Get differential expression results
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resdata$Gene <- as.character(resdata$Gene)
head(resdata)
hist(log10(res$pvalue), breaks=50, col="grey")
DESeq2::plotMA(dds)

# Volcano plot
## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.0005, topsig=10, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padjlfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
  with(subset(res, padjlfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  if (labelsig) {
    with(subset(res, -log10(pvalue) != Inf), text(log2FoldChange[1:topsig], -log10(pvalue[1:topsig]), labels=Gene, cex=textcx, col="blue"))
    with(subset(res, -log10(pvalue) != Inf), points(log2FoldChange[1:topsig], -log10(pvalue[1:topsig]), col="blue", pch=20, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both", "most significant"), pch=20, col=c("red","orange","green","blue"))
}
options(repr.plot.height=6, repr.plot.width=6)
volcanoplot(resdata, lfcthresh=2, sigthresh=0.00005, textcx=.8, topsig=5)
