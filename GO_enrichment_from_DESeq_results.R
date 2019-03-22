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
                         #FA_filename = "EveBCdTP_AnnotationTableFull.txt", #in case GO annotations aren't included
                          logFCthreshold = 3, #change to any FC threshold
                          padj.threshold = 0.001, #change to any adjusted p-value threshold
                          writeGenes = F) { #also write all genes for every term
  #read full data set
  full_list <- read.csv(filename, sep = sep)
  #rename variables we'll need later (probably...)
  #names(full_list)[3] <- "Sequence.Name"
  #names(full_list)[7] <- "Best diamond hit"
  #annotation <- read.delim(FA_filename) #don't need it if the file contains annotation
  #withAnnotation <- merge(x = annotation, y=full_list[,c(2,3,4,7)], all.x=F, all.y=F , by="Sequence.Name")
  #get de genes
  #either upregs
  if (upORdown == "up") {
    de <- full_list[full_list$log2FoldChange > logFCthreshold & full_list$padj < padj.threshold,]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "gene"] 
  }
  #or downregs
  if (upORdown == "down") {
    de <- full_list[full_list$log2FoldChange < -1 * logFCthreshold & full_list$padj < padj.threshold,]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "gene"]
  }  
  
  #get GO terms
  #if (gocat == "BP") BP <- withAnnotation[,c("Sequence.Name", "GO.Biological.Process")]
  #if (gocat == "MF") BP <- withAnnotation[,c("Sequence.Name", "GO.Molecular.Function")]
  if (gocat == "BP") BP <- full_list[,c("gene", "GO.Biological.Process")]
  if (gocat == "MF") BP <- full_list[,c("gene", "GO.Molecular.Function")]
  if (gocat == "CC") BP <- full_list[,c("gene", "GO.Cellular.Component")]
  
  #get only non-empty ones (I just cannot analyze all the rest)
  if (gocat == "BP") BPGO <- BP[BP$GO.Biological.Process != "-",]
  if (gocat == "MF")BPGO <- BP[BP$GO.Molecular.Function != "-",]
  if (gocat == "CC") BPGO <- BP[BP$GO.Cellular.Component != "-",] 
  #get all GO terms in machine-readable way (in a list + without names, only numbers)
  
  if (gocat == "BP") GOs <- strsplit(BPGO$GO.Biological.Process, split = "| ", fixed = T)
  if (gocat == "MF") GOs <- strsplit(BPGO$GO.Molecular.Function, split = "| ", fixed = T)
  if (gocat == "CC") GOs <- strsplit(BPGO$GO.Cellular.Component, split = "| ", fixed = T)
  names(GOs) <- BPGO$gene
  
  GOsTop <- lapply(X = GOs, function(x) gsub(" .*", "", x)) #remove human-readable name
  #get DE genes for the object
  DElist <- factor(as.integer(names(GOsTop) %in% de))
  
  ##exit if no DE genes
  if (length(levels(DElist)) < 2) {
    message("No DE genes to look at")
    return()}
  
  names(DElist) <- names(GOsTop)
  #construct a GOdat object (topGO package)
  GOdata <- new("topGOdata", ontology = gocat, allGenes = DElist, annot = annFUN.gene2GO, gene2GO = GOsTop)
  f <- runTest(GOdata, algorithm = "elim", statistic = "fisher") 
  #from the manual: We like to interpret the p-values returned by these methods as corrected or not affected by multiple testing.    
  signif_at_0.001 <- sum(f@score < 0.001)
  signif_at_0.05 <- sum(f@score < 0.05)
  allRes <- GenTable(object = GOdata, f, numChar = 100, topNodes = signif_at_0.05) #topNodes = signif_at_0.001 #nope! here without sorting
  allRes <- allRes[allRes$Significant > 1, ] #at least two genes. 
  
  if (writeGenes & nrow(allRes)) {
    allRes$Genes <- NA
    for (i in 1:length(allRes$Genes)) {
        temp <- genesInTerm(GOdata, allRes$GO.ID[i])[[1]]
        tempde <- temp[temp %in% de]
        namestempde <- full_list[full_list$gene %in% tempde, "best.nr.hit.diamond"]
        allRes$Genes[i] <- paste(namestempde, collapse = ", ")
    }
  }
  
  #output
  dir.name <- paste0("GO_FC", logFCthreshold, "_padj_", padj.threshold, "_", gocat)
  if (!dir.name %in% dir()) dir.create(dir.name)
  #full table
  names(allRes)[6] <- "p-value"
  write.csv(allRes, paste0(dir.name, "/", filename, "GO", upORdown, "_all", ".csv"))
  
  Res <- allRes[1:signif_at_0.001, ]  #allRes$`p-value` < 0.001, doesn't work properly
  
  write.csv(Res, paste0(dir.name, "/", filename, "GO", upORdown, ".csv"))
  return(allRes)
}

### Visualization (for particular case, but can be adjusted)
### A small helper function
process_GO_terms <- function(filename, species, condition, upORdown = "up") {
  if(file.exists(filename)) {
    table <- read.csv(filename)
    table5 <- table[table$Annotated > 2 & table$Significant > 1, c(2,3,7)] #at least 3 genes for a GO term # table$Annotated > 4
    #table5 <- table[table$Annotated > 2, c(2,3,7)] #at least 3 genes for a GO term # table$Annotated > 4
    #table5 <- table[, c(2,3,7)] #at least 3 genes for a GO term # table$Annotated > 4
    names(table5) <- c("GO.ID", "Term", paste0(species, condition))
    return(table5)
  }
  placeholder <- data.frame(GO.ID = NA, Term = NA, V3 = NA)
  names(placeholder)[3] <- paste0(species, condition)
  return(placeholder)
}

### A function for making GO term plots
gogo <- function(condition, upORdown = "up", suffix = "annot") {
  species = "Eve" #Eve
  filename <- paste0("Eve", condition, "03", "_", suffix, ".csvGO", upORdown, ".csv")
  eve3 <- process_GO_terms(filename, species, paste0(condition, "03"), upORdown = upORdown)
  filename <- paste0("Eve", condition, "24", "_", suffix, ".csvGO", upORdown, ".csv")
  eve24 <- process_GO_terms(filename, species, paste0(condition, "24"), upORdown = upORdown)
  species = "Ecy" #Ecy
  filename <- paste0("Ecy", condition, "03", "_", suffix, ".csvGO", upORdown, ".csv")
  ecy3 <- process_GO_terms(filename, species, paste0(condition, "03"),upORdown = upORdown)
  filename <- paste0("./Ecy", condition, "24", "_", suffix, ".csvGO", upORdown, ".csv")
  ecy24 <- process_GO_terms(filename, species, paste0(condition, "24"), upORdown = upORdown)
  species = "Gla" #Gla
  filename <- paste0("./Gla", condition, "03", "_", suffix, ".csvGO", upORdown, ".csv")
  gla3 <- process_GO_terms(filename, species, paste0(condition, "03"), upORdown = upORdown)
  filename <- paste0("./Gla", condition, "24", "_", suffix, ".csvGO", upORdown, ".csv")
  gla24 <- process_GO_terms(filename, species, paste0(condition, "24"), upORdown = upORdown)
  
  df <- join_all(list(eve3, eve24, ecy3, ecy24, gla3, gla24), by = c("GO.ID", "Term"), type = "full")
  
  
  df <- df[complete.cases(df$Term), ]
  
  df$Term <- ifelse(nchar(df$Term) > 33, 
                    paste0(substr(df$Term, 1, 30), "...", 
                           substr(df$Term, nchar(df$Term)-9, nchar(df$Term))),
                    df$Term)
  
  
  nnas <- apply(df, 1, function(x) sum(is.na(x)))
  names(nnas) <- 1:nrow(df)
  sortorder <- as.numeric(names(sort(nnas)))
  df <- df[sortorder, ]
  
  melted <- melt(df, id.vars = c("GO.ID", "Term"))
  names(melted) <- c("GO.ID", "Term", "Sample", "p.value")
  
  #cut toooo long names
  
  
  #make term factor
  melted$Term <- factor(melted$Term, levels = rev(df$Term))
  
  p <- melted$p.value
  p_value <- ifelse(p < 0.000001, "p < 10e-6", ifelse (p < 0.00001, "p < 10e-5", ifelse(p < 0.0001, "p < 10e-4", "p < 10e-3")))
  
  if (upORdown == "up") palette <- rev(brewer.pal(11, "PiYG")[1:4])
  if (upORdown == "down") palette <- brewer.pal(11, "PiYG")[8:11]
  
  p <- ggplot(melted, aes(x = Sample, y = Term, fill = p_value)) +
    geom_tile(aes(width=0.9, height=0.7), size=2) +
    scale_color_manual("black") +
    theme_classic() + 
    scale_fill_manual(values = palette) + 
    xlab("Species / treatment") + 
    scale_x_discrete(position = "top") + 
    geom_hline(yintercept = sum(nnas == 5) + 0.5, linetype = "dotted")
  
  width <- (max(nchar(df$Term)) + 72) / 13
  height = (length(df$Term) + 5) / 7
  
  ggsave(paste0("./", condition, "_", upORdown, ".png"), width = width, height = height)
  ggsave(paste0("./", condition, "_", upORdown, ".svg"), width = width, height = height)
  
  print(p)
  return(p)
}
### The same but now for time series
gogots <- function(upORdown = "up", species) {
  
  filename <- paste0("./", species, "T12", "_annot.csvGO", upORdown, ".csv")
  T12 <- process_GO_terms(filename, species, "T12", upORdown = upORdown)
  filename <- paste0("./", species, "T18", "_annot.csvGO", upORdown, ".csv")
  T18 <- process_GO_terms(filename, species, "T18", upORdown = upORdown)
  filename <- paste0("./", species, "T24", "_annot.csvGO", upORdown, ".csv")
  T24 <- process_GO_terms(filename, species, "T24", upORdown = upORdown)
  
  filename <- paste0("./", species, "LT1003", "_annot.csvGO", upORdown, ".csv")
  LT1003 <- process_GO_terms(filename, species, "LT1003", upORdown = upORdown)
  filename <- paste0("./", species, "LT1024", "_annot.csvGO", upORdown, ".csv")
  LT1024 <- process_GO_terms(filename, species, "LT1024", upORdown = upORdown)
  
  df <- join_all(list(T12, T18, T24, LT1003, LT1024), by = c("GO.ID", "Term"), type = "full") 
  df <- df[complete.cases(df$Term), ]
  
  nnas <- apply(df, 1, function(x) sum(is.na(x)))
  names(nnas) <- 1:nrow(df)
  sortorder <- as.numeric(names(sort(nnas)))
  df <- df[sortorder, ]
  
  melted <- melt(df, id.vars = c("GO.ID", "Term"))
  names(melted) <- c("GO.ID", "Term", "Sample", "p.value")
  
  melted$Term <- factor(melted$Term, levels = rev(df$Term))
  
  p <- melted$p.value
  p_value <- ifelse(p < 0.000001, "p < 10e-6", ifelse (p < 0.00001, "p < 10e-5", ifelse(p < 0.0001, "p < 10e-4", "p < 10e-3")))
  
  if (upORdown == "up") palette <- rev(brewer.pal(11, "PiYG")[1:4])
  if (upORdown == "down") palette <- brewer.pal(11, "PiYG")[8:11]
  
  p <- ggplot(melted, aes(x = Sample, y = Term, fill = p_value)) +
    geom_tile(aes(width=0.9, height=0.7), size=2) +
    scale_color_manual("black") +
    theme_classic() + 
    scale_fill_manual(values = palette) + 
    xlab("Species / treatment") + 
    scale_x_discrete(position = "top") + 
    geom_hline(yintercept = sum(nnas == 4) + 0.5, linetype = "dotted") +
    geom_vline(xintercept = 3.5, col = "grey")
  
  width <- (max(nchar(df$Term)) + 60) / 13
  height = (length(df$Term) + 5) / 7
  
  ggsave(paste0(species, "temp_gradual", "_", upORdown, ".png"), width = width, height = height)
  
  print(p)
  return(p)
}


###############################################################################
### Main part (example usage)
### GO term processing
### This is our main set
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
### Change here for another directory
#setwd("~/Documents/Paper1_stresses/data_tables/")
#setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/")
setwd("/run/media/polina/082A-8B6C/Paper1_stresses/data_tables/interspecies/")

setwd("~/Documents/Paper1_stresses/data_tables/decont_each_by_own/")

for (file in dir()) {
  if(grepl("csv", file)) {
    # ## Biological process, logFC cutoff of 3
    #  GOenrichment(filename = file, writeGenes = T)
    #  GOenrichment(filename = file, upORdown = "down", writeGenes = T)
    # # ## Molecular function, logFC cutoff of 3
    #  GOenrichment(filename = file, gocat = "MF", writeGenes = T)
    #  GOenrichment(filename = file, upORdown = "down", gocat = "MF", writeGenes = T)
    # ## Biological process, logFC 1
      GOenrichment(filename = file, gocat = "BP", logFCthreshold = 3, writeGenes = T)
      GOenrichment(filename = file, upORdown = "down", logFCthreshold = 3, writeGenes = T)
    # ## Molecular function, logFC 1
     #GOenrichment(filename = file, gocat = "MF", logFCthreshold = 3, writeGenes = T)
     #GOenrichment(filename = file, upORdown = "down", gocat = "MF", logFCthreshold = 3, writeGenes = T)
    ## Cellular component (just in case)
    #GOenrichment(filename = file, gocat = "CC", logFCthreshold = 3, writeGenes = T)
    #GOenrichment(filename = file, gocat = "CC", upORdown = "down", logFCthreshold = 3, writeGenes = T)
    
    
    #GOenrichment(filename = file, gocat = "BP", logFCthreshold = 1, writeGenes = T, padj.threshold = 0.05)
    #GOenrichment(filename = file, upORdown = "down", logFCthreshold = 1, writeGenes = T, padj.threshold = 0.05)
    # ## Molecular function, logFC 1
    GOenrichment(filename = file, gocat = "MF", logFCthreshold = 3, writeGenes = T)
    GOenrichment(filename = file, upORdown = "down", gocat = "MF", logFCthreshold = 3, writeGenes = T)
    ## Cellular component (just in case)
    GOenrichment(filename = file, gocat = "CC", logFCthreshold = 3, writeGenes = T)
    GOenrichment(filename = file, gocat = "CC", upORdown = "down", logFCthreshold = 3, writeGenes = T)
    
      }}




setwd("~/Documents/Paper2_time_series/compare_controls/data_tables/")
### Change here for another directory
setwd("~/Documents/Paper1_stresses/data_tables/interspecies/new")

setwd("~/Documents/Paper2_time_series/compare_controls/data_tables/Only_adjacent/")

for (file in dir()) {
  if(grepl("csv", file)) {
    GOenrichment(filename = file, gocat = "BP", logFCthreshold = 3, writeGenes = T)
    GOenrichment(filename = file, upORdown = "down", logFCthreshold = 3, writeGenes = T)

  }}




###############################################################################
##Let's now draw some illustrations##########
## manual switch here, which is bad. But easy and convenient
## Choose directory
# setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC1_padj_0.001_MF/")
# setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC1_padj_0.001_CC/")
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC3_padj_0.001_BP/")
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC3_padj_0.001_MF/")
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC3_padj_0.001_CC/")

  upp <- gogo("LT10", "up")
  downp <- gogo("LT10", "down")
  upcd <- gogo("Cd", "up")
  downcd <- gogo("Cd", "down")
  upac <- gogo("PB", "up")
  downac <- gogo("PB", "down")
  upph <- gogo("Ph", "up")
  downph <- gogo("Ph", "down")
  upacph <- gogo("PhB", "up")
  downacph <- gogo("PhB", "down")
  
setwd("~/Documents/Paper1_stresses/data_tables/decont_each_by_own/GO_FC3_padj_0.001_BP/")
upp <- gogo("LT10", "up", suffix = "decont")
downp <- gogo("LT10", "down", suffix = "decont")
upcd <- gogo("Cd", "up", suffix = "decont")
downcd <- gogo("Cd", "down", suffix = "decont")


### For full assemblies we need a switch, as the files are named differently
setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/GO1BP/")

 upp <- gogo("LT10", "up", suffix = "metazoa")
 downp <- gogo("LT10", "down", suffix = "metazoa")
 upcd <- gogo("Cd", "up", suffix = "metazoa")
 downcd <- gogo("Cd", "down", suffix = "metazoa")
 upac <- gogo("PB", "up", suffix = "metazoa")
 downac <- gogo("PB", "down", suffix = "metazoa")
 upph <- gogo("Ph", "up", suffix = "metazoa")
 downph <- gogo("Ph", "down", suffix = "metazoa")

###And time series (just to have the figures ready)
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO_FC1_padj_0.001_BP/")

gogots(species = "Eve")
gogots(species = "Ecy")
gogots(species = "Gla")

gogots(species = "Eve", upORdown = "down")
gogots(species = "Ecy", upORdown = "down")
gogots(species = "Gla", upORdown = "down")

setwd("/homes/brauerei/polina/Documents/Paper1_stresses/data_tables/interspecies")
GOenrichment("./EcyGlaB03_vs_B03_annot.csv", writeGenes = T)
GOenrichment("./EcyGlaB03_vs_B03_annot.csv", upORdown = "down", writeGenes = T)
GOenrichment("./GlaEcyB03_vs_B03_annot.csv", writeGenes = T)
GOenrichment("./GlaEcyB03_vs_B03_annot.csv", writeGenes = T, upORdown = "down")





################3
filename <- paste0("./", species, "T12", "_annot.csvGO", upORdown, ".csv")
T12 <- process_GO_terms(filename, species, "T12", upORdown = upORdown)
filename <- paste0("./", species, "T18", "_annot.csvGO", upORdown, ".csv")
T18 <- process_GO_terms(filename, species, "T18", upORdown = upORdown)
filename <- paste0("./", species, "T24", "_annot.csvGO", upORdown, ".csv")
T24 <- process_GO_terms(filename, species, "T24", upORdown = upORdown)

filename <- paste0("./", species, "LT1003", "_annot.csvGO", upORdown, ".csv")
LT1003 <- process_GO_terms(filename, species, "LT1003", upORdown = upORdown)
filename <- paste0("./", species, "LT1024", "_annot.csvGO", upORdown, ".csv")
LT1024 <- process_GO_terms(filename, species, "LT1024", upORdown = upORdown)

df <- join_all(list(T12, T18, T24, LT1003, LT1024), by = c("GO.ID", "Term"), type = "full") 
df <- df[complete.cases(df$Term), ]

nnas <- apply(df, 1, function(x) sum(is.na(x)))
names(nnas) <- 1:nrow(df)
sortorder <- as.numeric(names(sort(nnas)))
df <- df[sortorder, ]

melted <- melt(df, id.vars = c("GO.ID", "Term"))
names(melted) <- c("GO.ID", "Term", "Sample", "p.value")

melted$Term <- factor(melted$Term, levels = rev(df$Term))

p <- melted$p.value
p_value <- ifelse(p < 0.000001, "p < 10e-6", ifelse (p < 0.00001, "p < 10e-5", ifelse(p < 0.0001, "p < 10e-4", "p < 10e-3")))

if (upORdown == "up") palette <- rev(brewer.pal(11, "PiYG")[1:4])
if (upORdown == "down") palette <- brewer.pal(11, "PiYG")[8:11]

p <- ggplot(melted, aes(x = Sample, y = Term, fill = p_value)) +
  geom_tile(aes(width=0.9, height=0.7), size=2) +
  scale_color_manual("black") +
  theme_classic() + 
  scale_fill_manual(values = palette) + 
  xlab("Species / treatment") + 
  scale_x_discrete(position = "top") + 
  geom_hline(yintercept = sum(nnas == 4) + 0.5, linetype = "dotted") +
  geom_vline(xintercept = 3.5, col = "grey")

width <- (max(nchar(df$Term)) + 60) / 13
height = (length(df$Term) + 5) / 7

ggsave(paste0(species, "temp_gradual", "_", upORdown, ".png"), width = width, height = height)

print(p)
return(p)