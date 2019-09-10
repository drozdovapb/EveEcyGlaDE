#https://www.bioconductor.org/packages/devel/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html

#FAQ
#https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute") 
#biocLite("GO.db") 
####install.packages("qvalue")
#install.packages(c("dynamicTreeCut",,"flashClust","Hmisc","WGCNA"))

library(WGCNA)
options(stringsAsFactors = FALSE); enableWGCNAThreads()

#Loading the data; WGCNA requires genes be given in the columns

library(DESeq2)
library(tximport)
library(readr)
library(ggplot2)

#Important variables that should be located outside the presumable function

dir <- "/media/drozdovapb/big/Research/Projects/DE/salmon/eve_deep"; setwd(dir)

setwd(dir)
source("../de_functions_new+decontamination.R")

dir <- "/media/drozdovapb/ADATA NH13/drozdovapb/Research/Projects/DE/salmon/eve_deep"; setwd(dir)

samples <- read.csv("./samples.csv")
#create a vector of paths to files
files <- file.path(dir, "gc", samples$salmon.folder, "quant.sf")
names(files) <- paste0("sample", 1:6)
#check whether everything is fine and all the files exist
all(file.exists(files))
###choose samples
thispecies <- "Eve"
theseconditions <- c("B3", "B6", "B12", "B18", "B24")
to_use <- which(samples$species==thispecies & samples$condition %in% theseconditions)
txi <- tximport(files[to_use], type = "salmon", txOut=TRUE)


    sampleTable <- data.frame(condition=substr(samples$condition[to_use], 1, 1)) #use only the 1st letter
    #rownames(sampleTable) <- colnames(txi$counts) #it is NA
    rownames(sampleTable) <- samples$sample[to_use]
    colnames(txi$counts) <- samples$sample[to_use]
    
    #    Error in DESeqDataSet(se, design = design, ignoreRank) : 
    #    design has a single variable, with all samples having the same value.
    #use instead a design of '~ 1'. estimateSizeFactors, rlog and the VST can then be used
    
    #dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~1)
    
    keep <- rowSums(counts(dds)) >= 100 #about 4 per sample
    dds <- dds[keep,]
    #80 to 20 Mb It's great but how would I now match back? It would be quite complicated
    
    #differential expression! 
    dds <- DESeq(dds)
    #res <- results(dds)
    #res$log10padj <- -log10(res$padj)
    
    #vsd <- vst(dds, blind=FALSE)
    vsdd <- getVarianceStabilizedData(dds)
    plotPCA(vsd, intgroup=c("condition")) #well, there are some outliers...
    
    vsd.t <- t(vsd)
    


    
    
    