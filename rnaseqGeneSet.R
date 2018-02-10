# rnaseqGeneTest.R
# Date: 29/01/2018
# Author: Gabriel Forn-Cuní
# Description: The objective of this script is to perform gene set testing in 
# the pipeline following the rnaseq analysis. 
#
# Input: the .dgeResultsAnnot object from the pipeline
# Outputs: A folder with gene set testing for each comparison.

# From http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html:
# There are two kinds of tests: self-contained and competitive gene set tests. 
# Self-contained tests, which include the ROAST procedure, ask the question 
# "Are the genes in the set/pathway differentially expressed as a whole?" 
# Competitive gene set tests, like goana and camera ask the question whether the 
# differentially expressed genes tend to be over-represented in the gene set, 
# compared to all the other genes in the experiment. These different questions 
# use different statistical methodology.

# FIXME: Implement ROAST and CAMERA and GSEADB

# Edit the following variables.
# General
workingDirectory <- "C:/Users/gabri/Dropbox/Research/Leiden/projects/dram1/rnaseqDram"

outputRnks <- FALSE
runFgsea <- TRUE

# Do not edit below the line
#
# ---------------------------------------------------------------------------- #
#

# Load Utils
setwd(workingDirectory)
source("scripts/rnaseqUtils.R")

# load Data from the previous DGE Analysis
if(!file.exists('DGE/.dgeResultsAnnot')) {
  stop("Please run the DGE Analysis and Annotation prior the Enrichment.")
}
load("DGE/.dgeResultsAnnot")

library('fgsea')
library('ggplot2')

# Create the output directory
dir.create("geneSetAnalysis")
dir.create("geneSetAnalysis/rnks")
setwd(paste0(getwd(),"/geneSetAnalysis"))

for (x in 1:length(dgeResults$nonSig)) {
  comparison <- names(dgeResults$nonSig[x])
  
  for (y in 1:length(dgeResults$nonSig[[x]])) {
    method <- names(dgeResults$nonSig[[x]][y])
    print(paste(comparison, "-", method))
    #df <- data.frame(dgeResults$nonSig[[x]][[y]][,c("Homolog", "log2FC", "padj")])
    df <- data.frame(dgeResults$nonSig[[x]][[y]][,c("Zfin", "log2FC", "padj")])
    df <- df[order(df$padj),]
    df <- deleteDuplicatesDataFrame(df = df, col = "Zfin")
    #df <- deleteDuplicatesDataFrame(df = df, col = "Homolog")
    
    #df$rank <- -log10(df$padj)*df$log2FC
    df$rank <- -log10(df$padj)*sign(df$log2FC)
    
    rankedDF <- df[,c("Zfin", "rank")]
    rankedDF <- rankedDF[complete.cases(rankedDF),]
    #rankedDF$Symbol <- toupper(rankedDF$Symbol)

    if (outputRnks) {
      outputFileName <- paste0("rnks/", comparison, ".", method, ".geneSet.sign.zfin.rnk")
      
      write.table(rankedDF, file = outputFileName, sep = "\t", 
                  row.names=FALSE, col.names = FALSE, quote = FALSE)
    }
    if (runFgsea) {
      dir.create(comparison)
      ranks <- rankedDF
      ranks <- setNames(ranks$rank, ranks$Zfin)
      for (gmtFile in list.files("gmt")){
        print(paste("Analyzing subset", gmtFile))
        geneSet <- strsplit(gmtFile, "[.]")[[1]][1]
        pathways <- gmtPathways(paste0("gmt/",gmtFile))
        fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
        fgseaResSig <- subset(fgseaRes, padj <= 0.05)
        
        if (dim(fgseaResSig)[1] > 0) {
          
          # collapsedPathways <- collapsePathways(fgseaResSig[order(pval)][padj < 0.01], 
          #                                       pathways, ranks)
          # mainPathways <- fgseaResSig[pathway %in% collapsedPathways$mainPathways][
          #   order(-NES), pathway]

          topPathwaysUp <- fgseaResSig[ES > 0][head(order(padj), n=10)][order(-NES), pathway]
          topPathwaysDown <- fgseaResSig[ES < 0][head(order(padj), n=10)][order(NES), pathway]
          topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
          
          topFileName <- paste0(comparison, ".", method, 
                                ".", geneSet ,".topPlot.pdf")
          
          pdf(topFileName, 11.7, 8.3)
          
          print(plotGseaTable(pathways[topPathways], ranks, fgseaResSig, 
                              gseaParam = 0.2))
          dev.off()

          for (pathway in topPathways) {
            pathwayShort <- substr(pathway,0,20)
            pathwayPlot <- paste0(comparison, "/", comparison, ".", 
                                  method,  ".", geneSet , ".", pathwayShort, ".pdf")
            
            pdf(pathwayPlot, 11.7, 8.3)
            print(plotEnrichment(pathways[[pathway]], ranks) + labs(title=pathway))
            dev.off()
          }
          
        }
        
      }
    }   
  }
}


