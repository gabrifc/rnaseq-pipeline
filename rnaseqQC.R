# rnaseqQC.R
# Date: 13/12/2017
# Author: Gabriel Forn-Cuní
# Description: The objective of this script is to perform a Quality Check before
#   analyzing the differential gene expression on a RNAseq set.
#
# Input: A mapped reads count file, and a sample grouping file in csv format.
# Outputs:MDS plots using edgeR, and other useful charts to perform QC.

# Edit the name of the working directory and csv.
workingDirectory <- "C:/Users/gabri/Dropbox/Research/Leiden/projects/rnaseqDram"
inputFile <- "mappedReads.csv"
sampleGrouping <- "sampleGrouping.csv"

# Do not edit below the line
#
# ---------------------------------------------------------------------------- #
#

# Prepare libraries
if (!require("RColorBrewer")) {
  install.packages('RColorBrewer')
}
if (!require("gplots")) {
  install.packages('gplots')
}
if (!require("edgeR")) {
  source("http://www.bioconductor.org/biocLite.R") 
  biocLite("edgeR")
}

# Prepare style for graphs
myPalette <- brewer.pal(11,"RdYlBu")
myColors <- colorRampPalette(myPalette)(15)

# Get the files
setwd(workingDirectory)
readData <- read.csv(inputFile, header = TRUE, sep = ";", row.names = 1)
sampleGrouping <- read.csv(sampleGrouping, sep=";", row.names = 1)

# Create the output directory
dir.create("QC")
setwd(paste0(workingDirectory,"/QC"))

# CPM filtering (for edgeR)
countsPerMillion <- cpm(readData)
minLibrarySize <- min(colSums(readData))
CPMTreshold <- 10/(minLibrarySize/1e6)
readCountAboveTreshold <- countsPerMillion > CPMTreshold
readsToKeep <- rowSums(readCountAboveTreshold) >= 3
filteredReads <- readData[readsToKeep,]

# Prepare custom objects
dgeEdgeR <- DGEList(filteredReads)
dgeEdgeRNF <- DGEList(readData)
dgeEdgeRNormalized <- calcNormFactors(dgeEdgeR)

# Library Sizes Plots
par(mar=c(6,5,3,1))
pdf("librarySizes.pdf", 11.7, 8.3)
barplot(dgeEdgeRNF$samples$lib.size, 
        col = myColors,
        names = colnames(dgeEdgeR),
        las = 2,
        cex.names = 0.7,
        xlab = "",
        ylab = list("Size", font = 2),
        main = list("Library Size before filtering", font=4))
dev.off()

pdf("filteredLibrarySizes.pdf", 11.7, 8.3)
barplot(dgeEdgeR$samples$lib.size, 
                       col = myColors,
                       names = colnames(dgeEdgeR),
                       las = 2,
                       cex.names = 0.7,
                       xlab = list("Library", font = 2),
                       ylab = list("Size", font = 2),
                       main = list("Filtered Library Size", font=4))
dev.off()

# Boxplots distribution
logsCPM <- cpm(dgeEdgeR, log = TRUE)
logsCPMNF <- cpm(dgeEdgeRNF, log = TRUE)
logsCPMNR <- cpm(dgeEdgeRNormalized, log = TRUE)

pdf("filteredLibraryDistribution.pdf", 11.7, 8.3)
boxplot(logsCPM, 
        col = myColors,
        las = 2,
        cex.names = 0.7,
        xlab = "", 
        ylab = list("log2 CPM", font = 2), 
        main = list("Boxplot of filtered library CPM distribution", 
                    font = 4))
abline(h=median(logsCPM))
dev.off()

pdf("libraryDistribution.pdf", 11.7, 8.3)
boxplot(logsCPMNF, 
        col = myColors,
        las = 2,
        cex.names = 0.7,
        xlab = "", 
        ylab = list("log2 CPM", font = 2), 
        main = list("Boxplot of library CPM distribution before filtering", 
                    font = 4))
abline(h=median(logsCPMNF))
dev.off()

pdf("filteredLibraryDistributionNormalized.pdf", 11.7, 8.3)
boxplot(logsCPMNR, 
        col = myColors,
        las = 2,
        cex.names = 0.7,
        xlab = "", 
        ylab = list("log2 CPM", font = 2), 
        main = list("Boxplot of filtered library CPM distribution after TMM Normalization", 
                    font = 4))
abline(h=median(logsCPM))
dev.off()

# plot MDS
levelColors <- colorRampPalette(myPalette)(length(levels(sampleGrouping[,1])))
MDSColors <- levelColors[sampleGrouping[,1]]
if (length(sampleGrouping == 2)) {
  levelSymbols <- 14 + c(1:length(levels(sampleGrouping[,2])))
  MDSSymbols <- levelSymbols[sampleGrouping[,2]]
} else {
  MDSSymbols <- 16
}
pdf("filteredMDSSymbols.pdf", 11.7, 8.3)
plotMDS(dgeEdgeR, 
        col = MDSColors,
        pch = MDSSymbols,
        main = list("MDS of the filtered samples", font = 4))
legend("topright",
       legend = levels(sampleGrouping[,1]),
       col = levelColors,
       pch = 16)
legend("bottomright",
       legend = levels(sampleGrouping[,2]),
       pch = levelSymbols)
dev.off()
pdf("MDSSymbols.pdf", 11.7, 8.3)
plotMDS(dgeEdgeRNF,
        col = MDSColors,
        pch = MDSSymbols,
        main = list("MDS of the samples before filtering", font = 4))
legend("topright",
       legend = levels(sampleGrouping[,1]),
       col = levelColors,
       pch = 16)
legend("bottomright",
       legend = levels(sampleGrouping[,2]),
       pch = levelSymbols)
dev.off()

pdf("filteredMDS.pdf", 11.7, 8.3)
plotMDS(dgeEdgeR, 
        col = MDSColors,
        main = list("MDS of the filtered samples", font = 4))
dev.off()
pdf("filteredMDSNormalized.pdf", 11.7, 8.3)
dgeEdgeRNormalized <- calcNormFactors(dgeEdgeR)
plotMDS(dgeEdgeRNormalized, 
        col = MDSColors,
        main = list("MDS of the filtered samples", font = 4))
dev.off()
pdf("MDS.pdf", 11.7, 8.3)
plotMDS(dgeEdgeRNF,
        col = MDSColors,
        main = list("MDS of the samples before filtering", font = 4))
dev.off()

# Heatmap
varianceEdgeR <- apply(logsCPM, 1, var)
mostVariableGenesEdgeR <- names(sort(varianceEdgeR, decreasing=TRUE))[1:500]
mostVariableGeneCountsEdgeR <- logsCPM[mostVariableGenesEdgeR,]
pdf("filteredSamplesHeatmap.pdf", 11.7, 8.3)
heatmap.2(mostVariableGeneCountsEdgeR,
        col = myColors, 
        main = list("Top 500 most variable genes across filtered samples", 
                    font = 4),
        scale = "row",
        trace = "none",
        labRow = "",
        srtCol = 60,
        cexCol = 0.7,
        density.info = 'histogram',
        dendrogram = 'column',
        key = FALSE)
dev.off()

varianceEdgeRNF <- apply(logsCPMNF, 1, var)
mostVariableGenesEdgeRNF <- names(sort(varianceEdgeRNF, decreasing=TRUE))[1:500]
mostVariableGeneCountsEdgeRNF <- logsCPMNF[mostVariableGenesEdgeRNF,]
pdf("samplesHeatmap.pdf", 11.7, 8.3)
heatmap.2(mostVariableGeneCountsEdgeRNF,
          col = myColors, 
          main = list("Top 500 most variable genes across filtered samples", 
                      font = 4),
          scale = "row",
          trace = "none",
          labRow = "",
          srtCol = 60,
          cexCol = 0.7,
          density.info = 'histogram',
          dendrogram = 'column',
          key = FALSE)
dev.off()

varianceEdgeRNR <- apply(logsCPMNR, 1, var)
mostVariableGenesEdgeRNR <- names(sort(varianceEdgeRNR, decreasing=TRUE))[1:500]
mostVariableGeneCountsEdgeRNR <- logsCPMNR[mostVariableGenesEdgeRNR,]
pdf("filteredSamplesHeatmapNormalized.pdf", 11.7, 8.3)
heatmap.2(mostVariableGeneCountsEdgeRNR,
          col = myColors, 
          main = list("Top 500 most variable genes across TMM normalized filtered samples", 
                      font = 4),
          scale = "row",
          trace = "none",
          labRow = "",
          srtCol = 60,
          cexCol = 0.7,
          density.info = 'histogram',
          dendrogram = 'column',
          key = FALSE)
dev.off()