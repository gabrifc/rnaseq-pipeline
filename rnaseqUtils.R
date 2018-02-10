# rnaseqUtilities.R
# Date: 15/12/2017
# Author: Gabriel Forn-Cuní
# Description: The objective of this script is to centralize helper functions 
#   and variables that are used across diferent stages of the rnaseq analysis
#   to help with code reuse.
#
# Input: nothing.
# Output: nothing.
#
# FIXME: Conditionally load EdgeR and DESeq2 libraries 
# FIXME: Add descriptions for edgeR DGE calling.
# FIXME: Delete hardcoded data.
# FIXME: Fix different types of conversions to only one.
# FIXME: Add descriptions for everything.
#
# ---------------------------------------------------------------------------- #
#                             Libraries & Variables
# ---------------------------------------------------------------------------- #
#

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
if (!require("DESeq2")) {
  source("http://www.bioconductor.org/biocLite.R") 
  biocLite("DESeq2")
}  
if (!require("VennDiagram")) {
  install.packages('VennDiagram')
}
if (!require("goseq")) {
  source("http://www.bioconductor.org/biocLite.R") 
  biocLite("goseq")
}
if (!require("org.Dr.eg.db")) {
  source("http://www.bioconductor.org/biocLite.R") 
  biocLite("org.Dr.eg.db")
}

# Prepare style for graphs
myPalette <- brewer.pal(11,"RdYlBu")
myColors <- colorRampPalette(myPalette)(15)
vennColors <- c("#2980b9", "#16a085", "#27ae60")
# to prevent logs from venn.diagram()
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Function exportPDF.
# Input: a graph and a filename.
# Returns: nothing.
# Outputs: a pdf of the graph.

exportPDF <- function(graphToExport, fileName) {
  pdf(fileName, 11.7, 8.3)
  graphToExport
  dev.off()
}

# Function CPMfilter.
# Input: A read count matrix.
# Returns: A filtered count matrix, without read counts below 10/L, where L is 
# the smallest library size in millions.

CPMfilter <- function(readData) {
  return(readData[rowSums(cpm(readData) > 10/(min(colSums(readData))/1e6)) >= 2,])
}

# Function createVenn
# Inputs: a list of objects, title, subtitle and tiff filename.
# Returns: nothing. 

createVenn <- function(list, title, subtitle, fileName) {
  venn.diagram(list, 
               fill = vennColors, 
               alpha = 0.6, 
               filename = fileName, 
               lwd = 0.5,
               main = title,
               sub = subtitle,
               cat.default.pos = "text",
               main.fontface = 4,
               sub.fontface = 2,
               cat.fontface = 2,
               fontfamily = "sans",
               main.fontfamily = "sans",
               sub.fontfamily = "sans",
               cat.fontfamily = "sans",
               main.cex = 1.5,
               sub.cex = 1.3);
}

# Function collapseGrouping.
# Input: a sample Grouping of several columns.
# Returns: a one column sample grouping made by joining the rest of columns.

collapseGrouping <- function(sampleGrouping) {
  numberOfConditions <- length(sampleGrouping)
  nameOfConditions <- colnames(sampleGrouping)
  
  print(paste("Found", 
              length(sampleGrouping), 
              "conditions:", 
              paste(colnames(sampleGrouping), collapse = ", "),
              ". Collapsing for pairwise analysis."))
  
  sampleGrouping_args <- c(sampleGrouping, sep="")
  collapsedSampleGrouping <- as.data.frame(do.call(paste, sampleGrouping_args))
  rownames(collapsedSampleGrouping) <- colnames(readData)
  
  colnames(collapsedSampleGrouping) <- paste(colnames(sampleGrouping), 
                                             collapse = "")
  return(collapsedSampleGrouping)
}

# Function createDesignMatrix. 
# Input: sampleGrouping (1 column)
# Returns: the sample levels as a factor and a designMatrix.

createDesignMatrix <- function(sampleGrouping, sampleNames) {
  group <- as.factor(do.call(paste, c(sampleGrouping, sep="")))
  designMatrix <- model.matrix(~ 0 + group)
  colnames(designMatrix) <- levels(group)
  rownames(designMatrix) <- sampleNames
  returnList <- list("group" = group,
                     "designMatrix" = designMatrix)
  return(returnList)
}

# Function commonEdgeR.
# Input: a readData count file, and a sample grouping.
# Output: a DGEListObject with the samples grouped, TMM-Normalized, and with the
#   dispersion estimated.

commonEdgeRPrep <- function(readData, grouping) {
  y <- DGEList(readData)
  y$samples$group <- grouping$group
  y <- calcNormFactors(y)
  
  # Estimate the dispersion and export a graphic.
  y <- estimateDisp(y, grouping$designMatrix)
  exportPDF(plotBCV(y), "edgeR/dispersion.pdf")
  
  return(y)
}

# Function: createEdgeRContrasts.
# Input: a list of all the contrasts that will be evaluated.
# Return: a vector of the contrasts in edgeR format.

createEdgeRContrasts <- function(allComparisons) {
  contrastVector <- vector('character')
  for (i in 1:length(allComparisons)) {
    contrastVector[i] <- paste0(allComparisons[[i]][1], 
                                "-", 
                                allComparisons[[i]][2])
  }
  return(contrastVector)
}

# Function: dgeDESeq2.
# Inputs: an unfiltered read count matrix, the sample grouping and a list of all
#   contrasts to export.
# Returns: a list of lists with the significant genes for each comparison

dgeDESeq2 <- function(readData, sampleGrouping, allComparisons, significantList, nonSignificantList) {
  if (!all(rownames(sampleGrouping) %in% colnames(readData))) {
    message <- paste("The rownames of", 
                     sampleGroupingFile, 
                     "do not coincide with the columns of",
                     inputFile)
    stop(message)
  }
  if (!all(rownames(sampleGrouping) == colnames(readData))) {
    message <- paste("The order of the", 
                     sampleGroupingFile, 
                     "rows do not coincide with the columns of",
                     inputFile)
    stop(message)
  }
  
  # Prepare data for DESeq2
  design <- as.formula(paste("~", 
                             paste(colnames(sampleGrouping), 
                                   collapse=" + ")))
  
  dds <- DESeqDataSetFromMatrix(countData = readData,
                                colData = sampleGrouping,
                                design = design)
  
  # Plot a new PCA
  exportPDF(plotPCA(varianceStabilizingTransformation(dds), 
                    intgroup = colnames(sampleGrouping)), "DESeq2/PCA.pdf")
  
  # Start DESeq
  dds <- DESeq(dds)
  
  # Calculate foldchange for all of the contrasts
  for (i in 1:length(allComparisons)) {
    
    comparison <- c(colnames(sampleGrouping),
                    allComparisons[[i]][1],
                    allComparisons[[i]][2])
    
    ddsResults <- results(dds, 
                          alpha = 0.05, 
                          contrast = c(comparison[1], 
                                       comparison[2], 
                                       comparison[3]))
    
    ddsResultsOrdered <- ddsResults[order(ddsResults$padj),]
    ddsResultsSig <- subset(ddsResultsOrdered, padj < 0.05)
    
    # Output the data
    outputFileName <- paste0("DESeq2/pairwiseComparisons/",
                             allComparisons[[i]][1], 
                             "vs", 
                             allComparisons[[i]][2])
    
    exportPDF(DESeq2::plotMA(ddsResults, ylim = c(-2,2)),
              paste0(outputFileName,".MA.pdf"))
    
    write.csv(as.data.frame(ddsResultsSig), 
              file=paste0(outputFileName,".all.csv"))
    
    write.csv(as.data.frame(subset(ddsResultsSig, log2FoldChange > 0)), 
              file=paste0(outputFileName,".up.csv"))
    
    write.csv(as.data.frame(subset(ddsResultsSig, log2FoldChange < 0)), 
              file=paste0(outputFileName,".down.csv"))
    
    # Prepare the return
    
    # Check
    if (length(significantList) < i) {
      significantList[[i]] <- list()
    }
    
    sigDF <- as.data.frame(ddsResultsSig[, c("log2FoldChange", "pvalue", "padj")])
    significantList[[i]]$df$DESeq2 <- sigDF
    #significantList[[i]]$all$DESeq2 <- ddsResultsSig$log2FoldChange
    #names(significantList[[i]]$all$DESeq2) <- rownames(ddsResultsSig)
    significantList[[i]]$all$DESeq2 <- rownames(ddsResultsSig)
    significantList[[i]]$up$DESeq2 <- rownames(subset(ddsResultsSig, 
                                                      log2FoldChange > 0))
    significantList[[i]]$down$DESeq2 <- rownames(subset(ddsResultsSig, 
                                                        log2FoldChange < 0))
    
    if (length(nonSignificantList) < i) {
      nonSignificantList[[i]] <- list()
    }
    
    nonSigDF <- as.data.frame(ddsResults[, c("log2FoldChange", "pvalue", "padj")])
    nonSignificantList[[i]]$DESeq2 <- nonSigDF 
    #nonSignificantList[[i]]$DESeq2 <- ddsResults$log2FoldChange
    #names(nonSignificantList[[i]]$DESeq2) <- rownames(ddsResults)
  }
  
  output <- list(sig = significantList, nonSig = nonSignificantList)
  return(output)

}

# edgeRQL analysis

dgeEdgeRQL <- function(y, designMatrix, contrastVector, significantList, nonSignificantList) {
  
  fit <- glmQLFit(y, designMatrix, robust=TRUE)
  
  exportPDF(plotQLDisp(fit), "edgeR/fitQL.pdf")
  
  for (i in 1:length(contrastVector)) {
    
    # Make it global to work
    comparison <<- contrastVector[i]
    
    contrast <- makeContrasts(comparison, levels = designMatrix)
    test <- glmQLFTest(fit, contrast = contrast)
    
    edgerResults <- topTags(test, 
                           n = Inf, 
                           adjust.method = "BH", 
                           sort.by = "PValue", 
                           p.value = 0.05)
    
    edgerResults <- as.data.frame(edgerResults)
    
    allNS <- topTags(test, 
                       n = Inf, 
                       adjust.method = "BH", 
                       sort.by = "PValue", 
                       p.value = 1)
    
    allNS <- as.data.frame(allNS)
    
    # Output
    
    outputFileName <- paste0("edgeR/edgeRQL/",
                             strsplit(contrastVector[i], "-")[[1]][1], 
                             "vs", 
                             strsplit(contrastVector[i], "-")[[1]][2])
    
    exportPDF(plotMD(test), paste0(outputFileName,".MD.pdf"))
    
    
    write.csv(edgerResults, 
              file=paste0(outputFileName,".all.csv"))
    
    write.csv(subset(edgerResults, logFC > 0), 
              file=paste0(outputFileName,".up.csv"))
    
    write.csv(subset(edgerResults, logFC < 0), 
              file=paste0(outputFileName,".down.csv"))
    
    # Return
    
    # Check
    if (length(significantList) < i) {
      significantList[[i]] <- list()
    }
    
    sigDF <- edgerResults[,c("logFC", "PValue", "FDR")]
      
    significantList[[i]]$df$edgeRQL <- sigDF
    significantList[[i]]$all$edgeRQL <- rownames(edgerResults)
    #significantList[[i]]$all$edgeRQL <- significant$logFC
    #names(significantList[[i]]$all$edgeRQL) <- rownames(significant)
    significantList[[i]]$up$edgeRQL <- rownames(subset(edgerResults, 
                                                       logFC > 0))
    significantList[[i]]$down$edgeRQL <- rownames(subset(edgerResults, 
                                                         logFC < 0))
    
    if (length(nonSignificantList) < i) {
      nonSignificantList[[i]] <- list()
    }
    
    nonSigDF <- allNS[, c("logFC", "PValue", "FDR")]
    nonSignificantList[[i]]$edgeRQL <- nonSigDF 
    
    #nonSignificantList[[i]]$edgeRQL <- allNS$logFC
    #names(nonSignificantList[[i]]$edgeRQL) <- rownames(allNS)
  }
  
  output <- list(sig = significantList, nonSig = nonSignificantList)
  return(output)
}

dgeEdgeRLRT <- function(y, designMatrix, contrastVector, significantList, nonSignificantList) {
  fit <- glmFit(y, designMatrix, robust=TRUE)
  
  exportPDF(gof(fit, plot=T), "edgeR/fitLRT.pdf")
  
  for (i in 1:length(contrastVector)) {
    
    # Make it global to work
    comparison <<- contrastVector[i]
    contrast <- makeContrasts(comparison, levels = designMatrix)
    test <- glmLRT(fit, contrast = contrast)
    
    edgerResults <- topTags(test, 
                           n = Inf, 
                           adjust.method = "BH", 
                           sort.by = "PValue", 
                           p.value = 0.05)
    
    edgerResults <- as.data.frame(edgerResults)
    
    allNS <- topTags(test, 
                     n = Inf, 
                     adjust.method = "BH", 
                     sort.by = "PValue", 
                     p.value = 1)
    
    allNS <- as.data.frame(allNS)
    # Output
    
    outputFileName <- paste0("edgeR/edgeRLRT/",
                             strsplit(contrastVector[i], "-")[[1]][1], 
                             "vs", 
                             strsplit(contrastVector[i], "-")[[1]][2])
    
    exportPDF(plotMD(test), paste0(outputFileName,".MD.pdf"))
    
    write.csv(edgerResults, 
              file=paste0(outputFileName,".all.csv"))
    
    write.csv(subset(edgerResults, logFC > 0), 
              file=paste0(outputFileName,".up.csv"))
    
    write.csv(subset(edgerResults, logFC < 0), 
              file=paste0(outputFileName,".down.csv"))
    
    # Return
    
    # Check
    if (length(significantList) < i) {
      significantList[[i]] <- list()
    }
    
    sigDF <- edgerResults[,c("logFC", "PValue", "FDR")]
    
    significantList[[i]]$df$edgeRLRT <- sigDF
    significantList[[i]]$all$edgeRLRT <- rownames(edgerResults)
    
    #significantList[[i]]$all$edgeRLRT <- significant$logFC
    #names(significantList[[i]]$all$edgeRLRT) <- rownames(significant)
    significantList[[i]]$up$edgeRLRT <- rownames(subset(edgerResults, 
                                                        logFC > 0))
    significantList[[i]]$down$edgeRLRT <- rownames(subset(edgerResults, 
                                                          logFC < 0))
    
    if (length(nonSignificantList) < i) {
      nonSignificantList[[i]] <- list()
    }
    
    nonSigDF <- allNS[, c("logFC", "PValue", "FDR")]
    nonSignificantList[[i]]$edgeRLRT <- nonSigDF 

    #nonSignificantList[[i]]$edgeRLRT <- allNS$logFC
    #names(nonSignificantList[[i]]$edgeRLRT) <- rownames(allNS)
  }

  output <- list(sig = significantList, nonSig = nonSignificantList)
  return(output)
  
}

# ahDb <- query(ah, pattern = c("Danio Rerio", "GRCz10", "EnsDb")) to get the ID
# Input: the ID of the EnsDB object.
# Output: A vector with the gene Length Data.
downloadGeneLengthData <- function(id) {
  # Update the gene lengths
  library(ensembldb)
  library(AnnotationHub)
  ah <- AnnotationHub()
  ahDb <- query(ah, "EnsDb")
  EnsDb <- ahDb[[id]]
  return(lengthOf(EnsDb, of="gene"))
}

orgDrConversion <- function(listOfGenes) {
  library("org.Dr.eg.db")
  # Converting Ensembl gene ids to entrez gene ids for goana using org.Dr.eg.db.
  print("Converting Ensembl Gene IDs to Entrez Gene IDs")
  print("Conversion using org.Dr.eg.db")
  entrezGenesDB <- mapIds(org.Dr.eg.db,
                          keys=listOfGenes,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
  
  # Delete the Entrez NA genes (no data available)
  print("Deleting Repeated and Uncomplete Data.")
  
  # Delete all rows with duplicate EntrezGene IDs
  dup.idx = which(duplicated(entrezGenesDB))
  entrezGenesDB = entrezGenesDB[-dup.idx]
  
  print(paste("org.Dr.eg.db - Got", length(entrezGenesDB), "hits out of", length(listOfGenes)))
  
  return(entrezGenesDB)
} 

bioMartConversion <- function(listOfGenes, 
                              attributes = attributes,
                              filters = filter) {
  # Converting Ensembl gene ids to entrez gene ids for goana using Biomart
  library("biomaRt")
  print("Conversion using bioMaRt")
  zfishMart = useMart("ensembl", dataset="drerio_gene_ensembl")
  entrezGenesBM <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), 
                         filters = "ensembl_gene_id",
                         values = listOfGenes, mart = zfishMart)
  
  # Delete the Entrez NA genes (no data available)
  print("Deleting Repeated and Uncomplete Data.")
  entrezGenesBM <- entrezGenesBM[complete.cases(entrezGenesBM),]
  
  egbm <- entrezGenesBM[,2]
  names(egbm) <- entrezGenesBM[,1]
  
  # Delete all rows with duplicate EntrezGene IDs
  dup.idx <- which(duplicated(egbm))
  egbm <- egbm[-dup.idx]
  
  print(paste("bioMaRt - Got", length(egbm), "hits out of", length(listOfGenes)))
  return(egbm)
}

convertEnsembltoEntrez <- function(listOfGenes) {
  
  # TODO: For gage, first add a new row to allGenes that is entrezGeneID. 
  # Then, do the analysis for each comparison.
  
  entrezGenesDB <- orgDrConversion(listOfGenes)
  entrezGenesBM <- bioMartConversion(listOfGenes)
  if (length(entrezGenesDB) >= length(entrezGenesBM)) {
    print("Keeping ID conversion with org.Dr.eg.db")
    lostData <- length(listOfGenes) - length(entrezGenesDB)
    lostDataPercentage <- lostData/length(listOfGenes)*100
    print(paste0("We are losing info for ", lostData, " genes, or ", 
                 lostDataPercentage, "% of our analysis."))
    return(entrezGenesDB)
  } else {
    print("Keeping ID conversion with bioMaRt")
    lostData <- length(listOfGenes) - length(entrezGenesBM)
    lostDataPercentage <- lostData/length(listOfGenes)*100
    print(paste0("We are losing info for ", lostData, " genes, or ", 
                 lostDataPercentage, "% of our analysis."))
    return(entrezGenesBM)
  }
}

# Function deleDuplicatesDataFrame.
# Input: a dataframe and a column with diplucates.
# Output: the dataframe without the duplicates.
deleteDuplicatesDataFrame <- function(df, col) {
  dup.idx <- which(duplicated(df[col]))
  return(df[-dup.idx,])
}

bioMartConversion <- function(listOfGenes, 
                              attributes = attributes,
                              filters = filters) {
  # Converting Ensembl gene ids to entrez gene ids for goana using Biomart
  library("biomaRt")
  print("Annotating using BioMaRt")
  zfishMart = useMart("ensembl", dataset="drerio_gene_ensembl")
  geneAnnot <- getBM(attributes = attributes, 
                     filters = filters,
                     values = listOfGenes, mart = zfishMart)
  
  return(geneAnnot)
}

bioMaRtOrthologs <- function(listOfGenes, 
                             attributes = attributes,
                             filters = filters) {
  # Converting Ensembl gene ids to entrez gene ids for goana using Biomart
  library("biomaRt")
  humanDB = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  zfishMart = useMart("ensembl", dataset="drerio_gene_ensembl")
  print("Searching Orthologs with BioMaRt")
  orthologs <- getLDS(attributes = filters, filters = filters,
                      values = listOfGenes, mart = zfishMart, 
                      attributesL = attributes, martL = humanDB)
  
  return(orthologs)
}


# Function bioMaRtAnnot
# Input: list of Ensembl IDs
# Output: annotated dataframe
bioMaRtAnnot <- function(listOfIDS) {
  a <- bioMartConversion(listOfIDS, attributes = annotIDS, filters = filterID)
  
  b <- bioMaRtOrthologs(listOfIDS, attributes = "hgnc_symbol", filters = filterID)
  
  # Put together. 1st: list of Ids as df
  annot <- as.data.frame(listOfIDS)
  rownames(annot) <- listOfIDS
  colnames(annot) <- "ensembl_gene_id"
  
  # Clean the annotation. Ping the error source
  print("Deleting duplicated rows, mapping only the 1st result for each gene.")
  a <- deleteDuplicatesDataFrame(a, "ensembl_gene_id")
  b <- deleteDuplicatesDataFrame(b, "Gene.stable.ID")
  annot <- merge(annot, a, by = "ensembl_gene_id", all.x = TRUE)
  annot <- merge(annot, b, by.x = "ensembl_gene_id", by.y = "Gene.stable.ID", 
                 all.x = TRUE)
  return(annot)
}

# Function myGeneAnnot.
# Input: a list of IDs
# Return: an annotated datafram
myGeneAnnot <- function(listOfIDS) {
  library("mygene")
  annot <- queryMany(listOfIDS, 
                     scopes="ensembl.gene", 
                     fields=c("symbol", "entrezgene", "name"), 
                     species="zebrafish")
  annot <- as.data.frame(annotMyGene)
  annot <- deleteDuplicatesDataFrame(annot, "query")
  return(annot)
}

# Function txdbAnnot
# Input: a list of IDs
# Return: an annotated datafram
txdbAnnot <- function(listOfIDS) {
  library("org.Dr.eg.db")
  print("Annotation using org.Dr.eg.db")
  
  cols <- c("ENTREZID", "SYMBOL", "ZFIN", "GENENAME")
  txdbAnnot <- select(org.Dr.eg.db, 
                      keys=listOfIDS, 
                      columns=cols, 
                      keytype="ENSEMBL",
                      multiVals="first")
  txdbAnnot <- deleteDuplicatesDataFrame(txdbAnnot, "ENSEMBL")
  return(txdbAnnot)
}

# Function slimGOs.
# Input: a list of GO categories, a GOSlim database, and the Ontology to look.
# Output: a dataframe of the slimmed GO terms.
slimGOs <- function(listOfGOs, slimCollection, ont) {
  goCollection <- GOCollection(listOfGOs)
  goSlimResults <- goSlim(goCollection, slimCollection, ont)
  goSlimResults <- subset(goSlimResults, Count > 0)
  return(goSlimResults)
}

# Function runGoSeq.
# Input: list of genes, geneLengthData for bias, statistical method.
# Output: a dataframe with the statistically significant enriched GOs.

runGoSeq <- function(genes, geneLengthData, goSeqMethod) {
  # Analysis
  # Fitting the Probability Weighting Function (PWF) to check if DE genes are
  # biased by lenght
  pwf=nullp(genes, "danRer10", "ensGene", bias.data = geneLengthData)
  
  # For every case below, we can only select the GO that wer are interested using
  # test.cats=c("GO:MF"): GO.MF=goseq(pwf,"danRer10","ensGene", test.cats=c("GO:MF"))
  
  # Analysis
  if (goSeqMethod == 'Wallenius') {
    GO = goseq(pwf, "danRer10", "ensGene")
  } else if (goSeqMethod == 'RandomSampling') {
    GO = goseq(pwf, "danRer10", "ensGene", method="Sampling", repcnt=1000)
  } else if (goSeqMethod == 'NoBias') {
    GO = goseq(pwf, "danRer10", "ensGene", method="Hypergeometric")
  }
  
  # Output
  output <- list()
  goResults <- subset(GO, p.adjust(over_represented_pvalue, method="BH")<.05)
  rownames(goResults) <- goResults$category
  unEnrichedGo <- subset(GO, p.adjust(under_represented_pvalue, method="BH")<.05)
  rownames(unEnrichedGo) <- unEnrichedGo$category
  
  output$enrichedGo <- goResults
  output$unEnrichedGo <- unEnrichedGo
  return(output)
}

# Function deleteDuplicatesVector.
# Input: a vector with diplucates.
# Output: the vetor without the duplicates.
deleteDuplicatesVector <- function(vector) {
  vector <- as.vector(vector)
  dup.idx = which(duplicated(vector))
  if (length(dup.idx) > 0) {
    return(vector[-dup.idx])
  } else {
    return(vector)
  }
}


# Function reportLosses.
# Input: 2 vectors
# Ouptut: none. Prints the difference of length between the vectors.
reportLosses <- function(original, final) {
  numberFinal <- length(final)
  numberOriginal <- length(original)
  diff <- numberOriginal - numberFinal
  diffPerc <- format(round(diff/numberOriginal*100, 2), nsmall = 2)
  print(paste0("Losing info for ", diff, " genes, or ", diffPerc, "%."))
}


# Function renderPathways.
# Input: list of Pathways and gene expression
# Output: none. Creates images for every pathway in the list.
renderPathways <- function(listOfPathways, foldChangeData, suffix, keggDataFolder) {
  library("pathview")
  # Limit for log2FC colors
  limit = list(gene=1.5, cpd=1) 
  # Low expression color
  low = list(gene="deepskyblue3",cpd="deepskyblue3")
  # Mid Expression color
  mid = list(gene="gray",cpd="gray")
  # High Expression color
  high = list(gene="gold1",cpd="gold1")
  pv.out.list <- sapply(listOfPathways, function(pid) pathview(
    gene.data = foldChangeData, pathway.id = pid, species = "dre", 
    limit=limit, low=low, mid=mid, high=high, min.nnodes = 0, node.sum="mean", 
    same.layer = F, out.suffix = suffix, kegg.dir = keggDataFolder))
}