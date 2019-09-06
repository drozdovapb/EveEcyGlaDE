rm(list = ls())


{library(DESeq2)
  library(tximport)
  library(readr)
  library(ggplot2)
  library(pheatmap)}
#install.packages("ggrepel")
library(ggrepel) #it's for labels only


######CHOOSE WORKING DIRECTORY##################################################
#some hard-coded possible options below#
#dir <- "/media/drozdovapb/big/Research/Projects/DE/salmon/eve_deep"; setwd(dir)
#dir <- "/media/main/sandbox/drozdovapb/DE/salmon/eve_deep"; setwd(dir)
#dir <- "/media/drozdovapb/Elements/transcriptome/DE/4-salmon/eve_deep/"
dir <- "/run/media/polina/Elements/transcriptome/DE/9-reduced/ecy"; setwd(dir)
#read samples table
samples <- read.csv("./samples.csv")
#create a vector of paths to files
files <- file.path(dir, "gc", samples$salmon.folder, "quant.sf")
names(files) <- paste0("sample", 1:6)
#check whether everything is fine and all the files exist
all(file.exists(files))
#################################

#function to avoid running the same code    
import_count <- function(thispecies, theseconditions) { #now, get colors right
  palette <- c("Blues", "Greens", "Oranges")
  ##colors <- c("blue4", "green4", "orange4")
  ##colors <- c("#0072B2", "#009E73", "#D55E00") #should be the color blind friendly panel but....
  colors <- c("#00E2D2", "#007656", "#D55E00") #I see it. Coblis sees it.
  
  ###Transcripts need to be associated with gene IDs for gene-level summarization.
  ###should I be fine with transcript-level summarization?
  
  #txi <- tximport(files, type = "salmon", txOut=TRUE)
  #it would never work at my laptop
  
  #this version only worked for one species
  #tocompare <- which(samples$species==thispecies & samples$condition %in% theseconditions)
  #this should work for multiple species
  tocompare <- which(samples$species %in% thispecies & samples$condition %in% theseconditions)
  txi <- tximport(files[tocompare], type = "salmon", txOut=TRUE)
  
  sampleTable <- data.frame(condition=samples$condition[tocompare])
  #rownames(sampleTable) <- colnames(txi$counts) #it is NA
  rownames(sampleTable) <- samples$sample[tocompare]
  colnames(txi$counts) <- samples$sample[tocompare]
  
  dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
  
  return(dds)
}





# #First: different conditions, same species####
# #We can try 24h only, as the reaction was more pronounced? 
# #Eve, all conditions?
# thispecies <- c("Eve", "Ecy")
# theseconditions <- c("B03", "B24h", "LT1003", "LT1024", "Cd03", "Cd24")
# 
# evedata <- import_count(thispecies, theseconditions)
# 
# #DE
# evet3 <- read.csv("./EveLT1003_annot_de.csv")[,2]
# evet24 <- read.csv("./EveLT1024_annot_de.csv")[,2]
# evecd3 <- read.csv("./EveCd03_annot_de.csv")[,2]
# evecd24 <- read.csv("./EveCd24_annot_de.csv")[,2]
# ecyt3 <- read.csv("./EcyLT1003_annot_de.csv")[,2]
# ecyt24 <- read.csv("./EcyLT1024_annot_de.csv")[,2]
# ecycd3 <- read.csv("./EcyCd03_annot_de.csv")[,2]
# ecycd24 <- read.csv("./EcyCd24_annot_de.csv")[,2]
# 
# TCdDE <- unique(c(evet3, evet24, evecd3, evecd24, ecyt3, ecyt24, ecycd3, ecycd24)) #>500
# keep <- names(evedata) %in% TCdDE
# evedataTCdDE <- evedata[keep,] #6
# 
# #  
# #plotting function
# pca_plot_me_all <- function(temperdata, blind = TRUE, nsub = 10, label = F, 
#                             plottitle = deparse(substitute(temperdata)), fontsize = 11) {
#   #variance stabilizing transformation:
#   vsd <- vst(temperdata, blind=blind, nsub=nsub)
#   #if it doesn't work, try this instead:
#   # #https://support.bioconductor.org/p/62246/#62250
#   # dds <- temperdata[ rowSums(counts(temperdata)) > 5, ]
#   # cts <- counts(dds)
#   # geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
#   # dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
#   # vsd <- vst(dds, blind=blind, nsub=nsub)
#   
#   pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
#   pcaData$species <- substr(pcaData$name, 1, 3)
#   percentVar <- round(100 * attr(pcaData, "percentVar"))
#   p <- ggplot(pcaData, aes(PC1, PC2, col = species)) +
#     geom_point(size=3, stroke = 1.5, aes(shape=condition, fill=condition)) +
#     coord_fixed() + 
#     scale_shape_manual(values = c(21, 21, 25, 25, 23, 23)) +
#     scale_fill_manual(values = rep(c("gray60", "gray40"), 3)) +
#     scale_color_manual(values = c("#56B4E9", "#009E73")) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#     theme_light(base_size = fontsize) + ggtitle(plottitle)
#   if (label==T) p <- p + geom_label_repel(label = pcaData$name, label.padding = 1)
#   plot(p)
#   return(p)
# }
# 
# eveall <- pca_plot_me_all(evedataTCdDE, nsub = 10, label = F, plottitle = "Eve, all conditions")
# 
# #ggsave("~/Eve_all_conditions.png", eveall, width = 8, height = 6)




























#Eve, all conditions?
thispecies <- c("Eve", "Ecy")
theseconditions <- c("B24h", "LT1024", "Cd24")

evedata <- import_count(thispecies, theseconditions)

#DE
evet3 <- read.csv("./EveLT1003_annot_de.csv")[,2]
evet24 <- read.csv("./EveLT1024_annot_de.csv")[,2]
evecd3 <- read.csv("./EveCd03_annot_de.csv")[,2]
evecd24 <- read.csv("./EveCd24_annot_de.csv")[,2]
ecyt3 <- read.csv("./EcyLT1003_annot_de.csv")[,2]
ecyt24 <- read.csv("./EcyLT1024_annot_de.csv")[,2]
ecycd3 <- read.csv("./EcyCd03_annot_de.csv")[,2]
ecycd24 <- read.csv("./EcyCd24_annot_de.csv")[,2]

TCdDE <- unique(c(evet24, evecd24, ecyt24, ecycd24)) #~500
keep <- names(evedata) %in% TCdDE
evedataTCdDE <- evedata[keep,] #6

#  
#plotting function
pca_plot_me_all <- function(temperdata, blind = TRUE, nsub = 10, label = F, 
                            plottitle = deparse(substitute(temperdata)), fontsize = 11) {
  #variance stabilizing transformation:
  #vsd <- vst(temperdata, blind=blind, nsub=nsub)
  #if it doesn't work, try this instead:
  # #https://support.bioconductor.org/p/62246/#62250
   dds <- temperdata[ rowSums(counts(temperdata)) > 5, ]
   cts <- counts(dds)
   geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
   dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
   vsd <- vst(dds, blind=blind, nsub=nsub)
  
  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  pcaData$species <- substr(pcaData$name, 1, 3)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- ggplot(pcaData, aes(PC1, PC2, col = species)) +
    geom_point(size=3, stroke = 1.5, aes(shape=condition, fill=condition)) +
    coord_fixed() + 
    scale_shape_manual(values = c(21, 25, 23)) +
    scale_fill_manual(values = rep(c("gray60"), 3)) +
    scale_color_manual(values = c("#56B4E9", "#009E73")) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_light(base_size = fontsize) + ggtitle(plottitle)
  if (label==T) p <- p + geom_label_repel(label = pcaData$name, label.padding = 1)
  plot(p)
  return(p)
}

eveall <- pca_plot_me_all(evedataTCdDE, nsub = 100, label = F, fontsize = 16, plottitle = "")
ggsave("~/PCA_TCd_by_ecy.svg", width = 200, height = 150, units = "mm")
