# rnaseqDGE.R
# Date: 13/12/2017
# Author: Gabriel Forn-Cuní
# Description: The objective of this script is to analyze differential gene 
# expression on a RNAseq set with different methods and compare them.
#
# Input: A mapped reads count file, and a sample grouping file in csv format.
# Outputs: DGE called by the different algorithms chosen below.
#
# FIXME: CPM filter number of libraries chosen. 
# FIXME: Specify the reference conditions and normalize to them.
# FIXME: Specify if EdgeR should normalize by library size (TMM).
# FIXME: Refactor code for edgeR and significant genes exporting.
#

# Edit the following variables.
workingDirectory <- "C:/Users/gabri/Dropbox/Research/Leiden/projects/dram1/rnaseqDram"
inputFile <- "mappedReadsQC.csv"
sampleGroupingFile <- "sampleGroupingQC.csv"
#baseConditions <- c("WT", "PBS")

# Algorithms to use
collapseConditions <- TRUE
useCPMfiltering <- TRUE
useEdgeRQL <- TRUE
useEdgeRLRT <- TRUE
useDESeq2 <- TRUE
compareAlgorithms <- TRUE

# Do not edit below the line
#
# ---------------------------------------------------------------------------- #
#

# Load Utils
setwd(workingDirectory)
source("scripts/rnaseqUtils.R")

# Get the files
readData <- read.csv(inputFile, header = TRUE, sep = ";", row.names = 1)
sampleGrouping <- read.csv(sampleGroupingFile, sep=";", row.names = 1)

# Create the output directory
dir.create("DGE")
setwd(paste0(workingDirectory,"/DGE"))

# Start.
if (collapseConditions) {
  sampleGrouping <- collapseGrouping(sampleGrouping)
}

# Make contrasts
allComparisons <- combn(levels(sampleGrouping[,1]), 2, simplify = FALSE)
contrastVector <- createEdgeRContrasts(allComparisons)
#allGenes <- rownames(readData)
#save(allGenes, file=".allGenes")

# Create the dge lists
significant <- list()
nonSignificant <- list()

# DESeq2 Analysis
if (useDESeq2) {
  print("Running DESeq2")
  # Prepare the output folders
  dir.create("DESeq2")
  dir.create("DESeq2/pairwiseComparisons")
  # Analyze the data
  output <- dgeDESeq2(readData, sampleGrouping, allComparisons, 
                           significant, nonSignificant)
  significant <- output$sig
  nonSignificant <- output$nonSig
} 

# edgeR Analysis
if (useEdgeRLRT || useEdgeRQL) {
  # Create the output directory
  dir.create("edgeR")
  
  # Filter by CPM if specified
  if (useCPMfiltering) {
    print("Running CPM filter")
    readData <- CPMfilter(readData)
  }
  
  print("Starting edgeR")
  filteredGenes <- rownames(readData)
  
  # Common edgeR setup
  grouping <- createDesignMatrix(sampleGrouping, colnames(readData))

  # Fix for loop bugs (don't know why)
  designMatrix <- grouping$designMatrix
  
  y <- commonEdgeRPrep(readData, grouping)
  
  if (useEdgeRQL) {
    print("Running edgeRQL")
    dir.create("edgeR/edgeRQL")
    output <- dgeEdgeRQL(y, grouping$designMatrix, contrastVector, significant, nonSignificant)
    significant <- output$sig
    nonSignificant <- output$nonSig
  }
  if (useEdgeRLRT) {
    print("Running edgeRLRT")
    dir.create("edgeR/edgeRLRT")
    output <- dgeEdgeRLRT(y, grouping$designMatrix, contrastVector, significant, nonSignificant)
    significant <- output$sig
    nonSignificant <- output$nonSig
  }
}

# Save the significant list for later use.
names(significant) <- contrastVector
names(nonSignificant) <- contrastVector

dgeResults <- list(sig = significant, nonSig = nonSignificant)
#save(significant, file = ".significantGenes")
#save(nonSignificant, file = ".nonSignificantGenes")

save(dgeResults, file = ".dgeResults")

# Comparison
if (compareAlgorithms){
  print("Creating Venn Plots")
  dir.create("comparison")
  
  for (i in 1:length(significant)) {
    comparison <- names(significant[i])
    commonTitle <- paste("Common genes in", comparison)
    for (sub in significant[i]) {
      # All
      subtitle <- "All Significant Genes"
      fileName <- paste0("comparison/", comparison, ".all.tiff")
      createVenn(sub$all, commonTitle, subtitle, fileName)
      
      # Up
      subtitle <- "Significantly Upregulated Genes"
      fileName <- paste0("comparison/", comparison, ".up.tiff")
      createVenn(sub$up, commonTitle, subtitle, fileName)
      
      # Down
      subtitle <- "Significantly Downregulated Genes"
      fileName <- paste0("comparison/", comparison, ".down.tiff")
      createVenn(sub$down, commonTitle, subtitle, fileName)
    }
  }
}

print("Finished DGE analysis, choose method before enrichment.")