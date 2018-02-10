# rnaseqEnrichment.R
# Date: 15/12/2017
# Author: Gabriel Forn-Cuní
# Description: The objective of this script is to analyze DGE in a rnaseq 
#   dataset (use rnaseqDGE to choose one), annotate the dataset and then preform
#   enrichment analysis of Gene Ontology, Kegg Pathways and other dataset to 
#   discover the underlying altered biological functions. 
#
# Input: A mapped reads count file, and a sample grouping file in csv format.
# Outputs: A lot of things
#
# FIXME: write description.
# FIXME: Change all hardcorded Data about Zebrafish.
# FIXME: Functionize goana and move to rnaseqUtils.

# Edit the following variables.
# General
workingDirectory <- "C:/Users/gabri/Dropbox/Research/Leiden/projects/dram1/rnaseqDram"
printReportLosses <- TRUE
# Input
# Choose DGE algorithm chosen. Currently supports DESeq2, EdgeRQL, and EdgeRLRT
#DGEAlgorithm <- "DESeq2"

# Output
doGOAnalysis <- FALSE
doPathwayAnalysis <- TRUE

# GO
# GOAlgorithms implemented are: "goseq", "goana", and "clusterProfiler"
GOAlgoritm <- "goseq"
useGOSlim <- TRUE
# For GOSeq
# Method for calling significant results ("Wallenius, RandomSampling, 
# or NoBias (for gene length))
goSeqMethod <- "Wallenius"
# Download updated Data from Ensembl? For assemblies newer than DanRer6.
updateGeneLengthData <- TRUE

# Pathway
pathwayAlgorithm <- "kegga"

visualizePathways <- TRUE 

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

# Create the output directory
dir.create("Enrichment")
setwd(paste0(getwd(),"/Enrichment"))

# To minimize the amount of times we search for data, get info from  nonSig 
# which should include all genes. We only need one comparison, as all have the
# same genes. If CPMFilter was used, the edgeR nonSig are smaller than DESeq2,
# so first we try to get the DESeq2 data, and if it fails, get the edgeR.

for (i in dgeResults$nonSig) {
  # i = conditions
  if (!is.null(i[["DESeq2"]])) {
    listOfIDS <- i$DESeq2$Ensembl
  } else {
    # Get the first one, whatever that is
    listOfIDS <- i[[1]]$Ensembl
  }
}

if (doGOAnalysis) {
  if (GOAlgoritm == "goseq") {
    dir.create("goseq")
  } else if (GOAlgoritm == "goana") {
    dir.create("goana")
  } 
  # Download GOSlim data?
  if (useGOSlim) {
    print("Downloading GOSlim data")
    library("GSEABase") 
    slimURL <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
    slimCollection <- getOBOCollection(slimURL)
  } 
  # Update the gene Length Data?
  if (GOAlgoritm == "goseq" && GOMethod != "NoBias") {
    if (updateGeneLengthData) {
      # AH57746 | Ensembl 90 EnsDb for Danio Rerio
      geneLengthData <- downloadGeneLengthData('AH57746')
      genesWithLengthData <- listOfIDS[listOfIDS %in% names(geneLengthData)]
      geneLengthData <- geneLengthData[names(geneLengthData) %in% genesWithLengthData]
    } else {
      geneLengthData = NULL
    }
  }
}

if (doPathwayAnalysis) {
  if (pathwayAlgorithm == "kegga") {
    dir.create("kegga")
  }
  # else if (pathwayAlgorithm == "gage") {
  #   print(paste("Downloading KEGG Pathway list for zebrafish", species))
  #   kegg.path <- kegg.gsets(species = "dre",id.type = "kegg")
  #   keggPathways <- kegg.path$kg.sets
  # }
  if (visualizePathways) {
    dir.create("pathwayRenders")
    dir.create("pathwayRenders/keggData")
    keggDataFolder <- "pathwayRenders/keggData"
  }
}

# Start the analysis loop.
for (x in 1:length(dgeResults$sig)) {
  comparison <- names(dgeResults$sig[x])

  for (y in 1:length(dgeResults$sig[[x]]$df)) {
    method <- names(dgeResults$sig[[x]]$df[y])
    print(paste(comparison, "-", method))
    
    if (doGOAnalysis) {
      if (GOAlgoritm == "goseq") {
        # Analyse with goseq
        print("Using goseq")
        significantGenes <- dgeResults$sig[[x]]$df[[y]]$Ensembl
        genes <- as.integer(genesWithLengthData %in% significantGenes)
        names(genes) <- genesWithLengthData
        
        goSeqOutput <- runGoSeq(genes, geneLengthData, goSeqMethod)
        goResults <- goSeqOutput$enrichedGo
        unEnrichedGO <- goSeqOutput$unEnrichedGo
        
        write.csv2(as.data.frame(goResults), 
                   file=paste0("goseq/", comparison, ".", method, 
                               ".goseq.GO.csv"))

        write.csv2(as.data.frame(unEnrichedGO),
                   file=paste0("goseq/", comparison, ".", method, 
                               ".goseq.unEnrichedGO.csv"))
        
        if (useGOSlim) {
          listOfGOs <- rownames(goResults)
          for (ont in c("MF", "BP", "CC")) {
            goSlimResults <- slimGOs(listOfGOs, slimCollection, ont)
            fileName <- paste0("goseq/", comparison, ".", method, 
                               ".goseq.GOSlim.csv")
            write.csv2(goResults, file = fileName)
          }
        }
        
      } else if (GOAlgoritm == "goana") {
        print("Using goana")
        
        # Get Entrez ID of significant genes and all the dataset
        significantGenes <- na.omit(dgeResults$sig[[x]]$df[[y]]$Entrez)
        significantGenes <- deleteDuplicatesVector(significantGenes)
        
        universeGenes <- na.omit(dgeResults$nonSig[[x]][[y]]$Entrez)
        universeGenes <- deleteDuplicatesVector(universeGenes)
        
        if(printReportLosses) {
          print("Significant genes:")
          reportLosses(dgeResults$sig[[x]]$df[[y]]$Entrez, significantGenes)
          print("Whole Dataset:")
          reportLosses(dgeResults$nonSig[[x]][[y]]$Entrez, universeGenes)
        }
        
        goResults <- goana(significantGenes, universe = universeGenes, 
                           species = "Dr", prior.prob = NULL, covariate=NULL)
        
        goResults <- subset(goResults, p.adjust(P.DE, method="BH")<.05)
        
        fileName <- paste0("goana/", comparison, ".", method, ".goana.GO.csv")
        
        write.csv2(goResults, file = fileName)
        
        if (useGOSlim) {
          listOfGOs <- rownames(goResults)
          goSlimResults <- data.frame()
          for (ont in c("MF", "BP", "CC")) {
            ontSlimResults <- slimGOs(listOfGOs, slimCollection, ont)
            if (dim(ontSlimResults)[1] > 0) {
              ontSlimResults$ont <- ont
            }
            goSlimResults <- rbind(goSlimResults, ontSlimResults)
            fileName <- paste0("goana/", comparison, ".", method, 
                               ".goana.GOSlim.csv")
            write.csv2(goSlimResults, file = fileName)
          }
        }
      } else if (GOAlgoritm == "clusterProfiler") {
        
      }
    } 
    if (doPathwayAnalysis) {
      if (pathwayAlgorithm == "kegga") {
        print("Analysing pathways with kegga.")
        # Get Entrez ID of significant genes and all the dataset
        significantGenes <- na.omit(dgeResults$sig[[x]]$df[[y]]$Entrez)
        significantGenes <- deleteDuplicatesVector(significantGenes)
        
        universeGenes <- na.omit(dgeResults$nonSig[[x]][[y]]$Entrez)
        universeGenes <- deleteDuplicatesVector(universeGenes)
        
        if(printReportLosses) {
          print("Significant genes:")
          reportLosses(dgeResults$sig[[x]]$df[[y]]$Entrez, significantGenes)
          print("Whole Dataset:")
          reportLosses(dgeResults$nonSig[[x]][[y]]$Entrez, universeGenes)
        }
        
        pathwayResults <- kegga(significantGenes, universe = universeGenes, 
                           species = "Dr")
        
        pathwayResults <- subset(pathwayResults, p.adjust(P.DE, method="BH")<.05)
        
        fileName <- paste0("kegga/", comparison, ".", method, ".kegga.csv")

        write.csv2(pathwayResults, file = fileName)
        
        if (visualizePathways) {
          # Visualization variables 
          
          foldChangeData <- dgeResults$nonSig[[x]][[y]][,c("Entrez", "log2FC")]
          foldChangeData <- deleteDuplicatesDataFrame(foldChangeData, "Entrez")
          foldChangeList <- foldChangeData[,2]
          names(foldChangeList) <- foldChangeData$Entrez
          
          # pathwayResults, extract the number of the pathway. e.g. from 
          # "path:dre00022" to 00022
          listOfPathways <- c(substr(c(rownames(pathwayResults)), 9, 13))
          
          suffix <- paste0(comparison,".",method)
          
          renderPathways(listOfPathways, foldChangeData = foldChangeList, 
                         suffix = suffix, keggDataFolder = keggDataFolder)
        }
        
      } 
      # else if (pathwayAlgorithm == "gage") {
      #   print("Analysing pathways with gage.")
      # }
    }
    
  }
}
