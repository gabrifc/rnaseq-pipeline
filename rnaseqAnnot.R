# rnaseqAnnotation.R
# Date: 24/01/2018
# Author: Gabriel Forn-Cuní
# Description: The objective of this annotate a previous zebrafish DGE Analysis 
# with things that we could need later: Entrez ID, Gene Symbol, ZFIN, full name,
# first Homolog human gene name.
#
# Input: .dgeResults list created by rnaseqDGE.R
# Outputs: The DGE lists annotated
#
# FIXME: Add gProfiler (https://biit.cs.ut.ee/gprofiler/page.cgi?apis) 
# FIXME: Add biodbnet.
# FIXME: Fix filling NAs.

# Edit the following variables.
# General
workingDirectory <- "C:/Users/gabri/Dropbox/Research/Leiden/projects/dram1/rnaseqDram"

annotIDS <- c("ensembl_gene_id", "entrezgene", "external_gene_name", 
              "zfin_id_id", "description")
filterID <- "ensembl_gene_id"

# Select method to convert ids from dbs: biomart, mygene, org.Dr.eg.db or all
conversionMethod <- "txdb"

# Do not edit below the line
#
# ---------------------------------------------------------------------------- #
#

# Load Utils
setwd(workingDirectory)
source("scripts/rnaseqUtils.R")

setwd(paste0(getwd(),"/DGE"))

# load Data from the previous DGE Analysis
if(!file.exists('.dgeResults')) {
  stop("Please run the DGE Analysis prior the Enrichment.")
}
load(".dgeResults")

dir.create("annotation")
setwd(paste0(getwd(),"/annotation"))

# To minimize the amount of times we search for data, convert the nonSig list,
# which should include all genes. We only need one comparison, as all have the
# same genes. If CPMFilter was used, the edgeR nonSig are smaller than DESeq2,
# so first we try to get the DESeq2 data, and if it fails, get the edgeR.

for (i in dgeResults$nonSig) {
  # i = conditions
  if (!is.null(i[["DESeq2"]])) {
    listOfIDS <- rownames(i$DESeq2)
  } else {
    # Get the first one, whatever that is
    listOfIDS <- rownames(i[[1]])
  }
}

# After testing multiple ways of annotating the data, I have come to believe 
# that txdb is the best method, in terms of speed and accuracy. I will therefore
# use txdb to annotate the majority of the data, and fill the missing data with
# biomart. 

print("Annotating using TxDB. Note that homologs are not yet supported.")
annot <- txdbAnnot(listOfIDS)
colnames(annot) <- c("Ensembl", "Entrez", "Symbol", "Zfin", "Name")

print("Missing the missing data with biomart")
annotBiomart <- bioMaRtAnnot(listOfIDS)
# Hack, move columns
colnames(annotBiomart) <- c("Ensembl", "Entrez", "Symbol", "Name", "Zfin",
                            "Homolog")
#annotBiomart <- annotBiomart[,colnames(annot)]

# Merge homolog data
homologs <- annotBiomart[,c("Ensembl", "Homolog")]
annot <- merge(annot, homologs, by = "Ensembl", all.x = TRUE)

# # Fill missing data in annot
# for(i in 1:dim(annot)[1]){
#   annot[i,is.na(annot[i,])] <- annotBiomart[i,]
# }

 
# Now we have a data frame with all the values that we will need, let's annotate
# the result dataframes.
for (x in 1:length(dgeResults$sig)) {
  for (y in 1:length(dgeResults$sig[[x]]$df)) {
    dgeResults$sig[[x]]$df[[y]] <- merge(dgeResults$sig[[x]]$df[[y]], annot, 
                                         by.x = "row.names", by.y = "Ensembl", 
                                         all.x = TRUE)
    colnames(dgeResults$sig[[x]]$df[[y]]) <- c("Ensembl", "log2FC", "pvalue", 
                                               "padj", "Entrez", "Symbol", 
                                               "Zfin", "Name", "Homolog")
  }
}

for (x in 1:length(dgeResults$nonSig)) {
  comparison <- names(dgeResults$nonSig[x])
  for (y in 1:length(dgeResults$nonSig[[x]])) {
    dgeResults$nonSig[[x]][[y]] <- merge(dgeResults$nonSig[[x]][[y]], annot, 
                                         by.x = "row.names", by.y = "Ensembl", 
                                         all.x = TRUE)
    colnames(dgeResults$nonSig[[x]][[y]]) <- c("Ensembl", "log2FC", "pvalue", 
                                               "padj", "Entrez", "Symbol", 
                                               "Zfin", "Name", "Homolog")
    method <- names(dgeResults$nonSig[[x]][y])
    fileName <- paste0(comparison,".",method,".all.csv")
    write.csv2(dgeResults$nonSig[[x]][[y]], file = fileName)
    fileName <- paste0(comparison,".",method,".sig.csv")
    write.csv2(dgeResults$sig[[x]]$df[[y]], file = fileName)
    
  }
}

setwd(paste0(workingDirectory,"/DGE"))

save(dgeResults, file = ".dgeResultsAnnot")

# # We have the complete list of IDs now, let's get a table for conversion.
# if (conversionMethod %in% c("biomart", "all")) {
#   annotBiomart <- bioMaRtAnnot(listOfIDS)
# } else if (conversionMethod %in% c("mygene", "all")) {
#   print("Annotating using MyGene. Note that homologs and ZFIN are not yet supported.")
#   annotMyGene <- myGeneAnnot(listOfIDS)
# } else if (conversionMethod %in% c("txdb", "all")) {
#   print("Annotating using TxDB Note that homologs are not yet supported.")
#   annotTxdb <- txdbAnnot(listOfIDS)
# }

# dir.create("annotation")
# setwd(paste0(getwd(),"/annotation"))
# 
# for (x in 1:length(dgeResults$nonSig)) {
#   for (y in 1:length(dgeResults$nonSig[[x]])) {
#     method <- names(dgeResults$nonSig[[x]][y])
#     fileName <- paste0(comparison,".",method,".all.csv")
#     write.csv2(dgeResults$nonSig[[x]][[y]], file = fileName)
#     fileName <- paste0(comparison,".",method,".sig.csv")
#     write.csv2(dgeResults$sig[[x]]$df[[y]], file = fileName)
#   }
# }
