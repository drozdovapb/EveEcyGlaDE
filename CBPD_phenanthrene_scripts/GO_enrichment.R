####This script serves to get overrepresented GO terms from DE (DESeq2) results

###reading options
options(stringsAsFactors = F)
### graphical options
#options(device='x11')
#grDevices::X11.options(type='cairo')
options(bitmapType='cairo')

###libraries
require(plyr)
require(reshape2)
require(ggplot2)
library(RColorBrewer)
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library(topGO)

###just in case: color palette
#scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73",
#"#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
#black orange sky_blue green yellow blue vermillion reddish_purple


###############################################################################
###Functions
###This function is for list of genes => enriched GO terms. Supports custom FC cutoff.
##testing set
#filename = "EcyT12_annot.csv"; sep = ","; upORdown = "up"; gocat = "BP"
#logFCthreshold = 1; padj.threshold = 0.001; writeGenes = T
#fa_filename = "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt"

GOenrichment <- function(filename, sep = ",", upORdown = "up", gocat = "BP", 
                         DEfilename = "", 
                         FA_filename = "", #in case GO annotations aren't included
                         logFCthreshold = 1, #change to any FC threshold
                         padj.threshold = 0.05, #change to any adjusted p-value threshold
                         writeGenes = T) { #also write all genes for every term
  full_list <- read.csv(filename, sep = sep)
  #rename variables we'll need later (probably...)
  full_list$"Sequence.Name" <- row.names(full_list)
  annotation <- read.delim(FA_filename)
  #merge tables
  withAnnotation <- merge(x = annotation, y=full_list, all.x=F, all.y=F , by="Sequence.Name")
  #get de genes
  defile <- read.csv(DEfilename)
  #either upregs
  if (upORdown == "up") {
    de <- defile[defile$salmon.DESeq2.FC > 0, "Gene"]
    ##get rid of NAs & get names only
    #de <- de[complete.cases(de$padj), "Sequence.Name"] 
  }
  #or downregs
  if (upORdown == "down") {
    de <- defile[defile$salmon.DESeq2.FC < 0, "Gene"]
  }  
  #get GO terms
  if (gocat == "BP") BP <- withAnnotation[,c("Sequence.Name", "GO.Biological.Process")]
  if (gocat == "MF") BP <- withAnnotation[,c("Sequence.Name", "GO.Molecular.Function")]
  if (gocat == "CC") BP <- withAnnotation[,c("Sequence.Name", "GO.Cellular.Component")]
  
  
  #get only non-empty ones (I just cannot analyze all the rest)
  if (gocat == "BP") BPGO <- BP[BP$GO.Biological.Process != "-",]
  if (gocat == "MF")BPGO <- BP[BP$GO.Molecular.Function != "-",]
  if (gocat == "CC")BPGO <- BP[BP$GO.Cellular.Component != "-",]
  #get all GO terms in machine-readable way (in a list + without names, only numbers)
  
  if (gocat == "BP") GOs <- strsplit(BPGO$GO.Biological.Process, split = "| ", fixed = T)
  if (gocat == "MF") GOs <- strsplit(BPGO$GO.Molecular.Function, split = "| ", fixed = T)
  if (gocat == "CC") GOs <- strsplit(BPGO$GO.Cellular.Component, split = "| ", fixed = T)
  names(GOs) <- BPGO$Sequence.Name
  
  GOsTop <- lapply(X = GOs, function(x) gsub(" .*", "", x)) #remove human-readable name
  #get DE genes for the object
  DElist <- factor(as.integer(names(GOsTop) %in% de))
  names(DElist) <- names(GOsTop)
  #construct a GOdat object (topGO package)
  GOdata <- new("topGOdata", ontology = gocat, allGenes = DElist, annot = annFUN.gene2GO, gene2GO = GOsTop)
  f <- runTest(GOdata, algorithm = "elim", statistic = "fisher") 
  #from the manual: We like to interpret the p-values returned by these methods as corrected or not affected by multiple testing.    
  signif_at_0.01 <- sum(f@score < 0.01)
  allRes <- GenTable(object = GOdata, f, topNodes = signif_at_0.01, numChar = 100)        
  
  if (writeGenes & nrow(allRes)) {
    allRes$Genes <- NA
    for (i in 1:length(allRes$Genes)) {
      temp <- genesInTerm(GOdata, allRes$GO.ID[i])[[1]]
      tempde <- temp[temp %in% de]
      namestempde <- withAnnotation[withAnnotation$Sequence.Name %in% tempde, "best.hit.to.nr"]
      allRes$Genes[i] <- paste(namestempde, collapse = ", ")
    }
  }
  
  #output
  dir.name <- paste0("GO_FC", logFCthreshold, "_padj_", padj.threshold, "_", gocat)
  if (!dir.name %in% dir()) dir.create(dir.name)
  #full table
  names(allRes)[6] <- "p-value"
  write.csv(allRes, paste0(filename, "GO", gocat, upORdown, "_all", ".csv"))
  ## i tried with subdirectories, but they do not work if the original file was not in the same folder
  signif_at_0.001 <- sum(f@score < 0.001)
  Res <- allRes[1:signif_at_0.001, ]  #allRes$`p-value` < 0.001, doesn't work properly
  
  write.csv(Res, paste0(filename, "GO", upORdown, ".csv"))
  return(allRes)
}

setwd("/home/drozdovapb/Research/Projects/DE/texts/Paper1_stresses/acetone_phenanthrene_story/CBPD_submit/1-R1/new_DE_comparison/")

#file1 <- "PB03_vs_B03.Ecycommon_DE.csv"
#for (file in dir()) {
#  if(grepl("csv", file)) {
#    # ## Biological process, logFC cutoff of 3
#    #  GOenrichment(filename = file, writeGenes = T)
#    #  GOenrichment(filename = file, upORdown = "down", writeGenes = T)
#    # # ## Molecular function, logFC cutoff of 3
#    #  GOenrichment(filename = file, gocat = "MF", writeGenes = T)
#    #  GOenrichment(filename = file, upORdown = "down", gocat = "MF", writeGenes = T)
#    # ## Biological process, logFC 1
#    GOenrichment(filename = file1, gocat = "BP", logFCthreshold = 1, writeGenes = T)
#  }}    
    

#filename <- "../new_DE_comparison/salmon/Ecy.isoform.counts.matrix.B24h_vs_PB24.DESeq2.DE_results"
#sep = "\t"
FA_filename = "~/Research/Projects/DE/annotation/EcyBCdTP1_cor_AnnotationTable.txt"
##Ecy, 3h, solvent
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", writeGenes = T,
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Ecycommon_DE.csv")
##Ecy, 3h, phenanthrene
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Ecycommon_DE.csv")

##Ecy, 24h, solvent
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Ecycommon_DE.csv")
##Ecy, 24h, phenanthrene
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Ecycommon_DE.csv")



FA_filename = "~/Research/Projects/DE/annotation/EveBCdTP1_cor_AnnotationTable.txt"
##Eve, 3h, solvent
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Evecommon_DE.csv")
##Eve, 3h, phenanthrene
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Evecommon_DE.csv")

##Eve, 24h, solvent
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Evecommon_DE.csv")
##Eve, 24h, phenanthrene
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Evecommon_DE.csv")





FA_filename = "~/Research/Projects/DE/annotation/GlaBCdTP1_cor_AnnotationTable.txt"
##Gla, 3h, solvent
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Glacommon_DE.csv")
##Gla, 3h, phenanthrene
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Glacommon_DE.csv")

##Gla, 24h, solvent
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Glacommon_DE.csv")
##Gla, 24h, phenanthrene
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "BP", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Glacommon_DE.csv")


##MF
################



##CC
################
##CC
FA_filename = "~/Research/Projects/DE/annotation/EcyBCdTP1_cor_AnnotationTable.txt"
##Ecy, 3h, solvent
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", writeGenes = T,
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Ecycommon_DE.csv")
##Ecy, 3h, phenanthrene
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Ecycommon_DE.csv")

##Ecy, 24h, solvent
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Ecycommon_DE.csv")
##Ecy, 24h, phenanthrene
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Ecycommon_DE.csv")
GOenrichment(filename = "salmon/Ecy.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Ecycommon_DE.csv")



FA_filename = "~/Research/Projects/DE/annotation/EveBCdTP1_cor_AnnotationTable.txt"
##Eve, 3h, solvent
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Evecommon_DE.csv")
##Eve, 3h, phenanthrene
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Evecommon_DE.csv")

##Eve, 24h, solvent
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Evecommon_DE.csv")
##Eve, 24h, phenanthrene
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Evecommon_DE.csv")
GOenrichment(filename = "salmon/Eve.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Evecommon_DE.csv")





FA_filename = "~/Research/Projects/DE/annotation/GlaBCdTP1_cor_AnnotationTable.txt"
##Gla, 3h, solvent
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB03_vs_B03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB03_vs_B03.Glacommon_DE.csv")
##Gla, 3h, phenanthrene
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph03_vs_PB03.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph03_vs_PB03.Glacommon_DE.csv")

##Gla, 24h, solvent
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "up", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.PB24_vs_B24h.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "PB24_vs_B24h.Glacommon_DE.csv")
##Gla, 24h, phenanthrene
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Glacommon_DE.csv")
GOenrichment(filename = "salmon/Gla.isoform.counts.matrix.Ph24_vs_PB24.DESeq2.DE_results", 
             upORdown = "down", gocat = "CC", sep = "\t", 
             FA_filename = FA_filename, DEfilename = "Ph24_vs_PB24.Glacommon_DE.csv")