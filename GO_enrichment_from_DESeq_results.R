####This script is for GO terms

options(stringsAsFactors = F)
library(topGO)

GOenrichment <- function(filename, sep = ",", upORdown = "up", gocat = "BP", 
                         FA_filename = "/homes/brauerei/polina/Documents/transcriptome_annotation/EveBCdTP_AnnotationTableFull.txt") {
  #read full data set
  full_list <- read.csv(filename, sep = sep)
  #rename variables we'll need later (probably...)
  #names(full_list)[3] <- "Sequence.Name"
  #names(full_list)[7] <- "Best diamond hit"
  #annotation <- read.delim(FA_filename) #don't need it if the file contains annotation
  ##merge tables
  #withAnnotation <- merge(x = annotation, y=full_list[,c(2,3,4,7)], all.x=F, all.y=F , by="Sequence.Name")
  #get de genes
  #either upregs
  if (upORdown == "up") {
    de <- full_list[full_list$log2FoldChange > 3 & full_list$padj < 0.001,]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "gene"] 
  }
  #or downregs
  if (upORdown == "down") {
    de <- full_list[full_list$log2FoldChange < -3 & full_list$padj < 0.001,]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "gene"]
  }  
  
  
  
  #get GO terms
  #if (gocat == "BP") BP <- withAnnotation[,c("Sequence.Name", "GO.Biological.Process")]
  #if (gocat == "MF") BP <- withAnnotation[,c("Sequence.Name", "GO.Molecular.Function")]
  if (gocat == "BP") BP <- full_list[,c("gene", "GO.Biological.Process")]
  if (gocat == "MF") BP <- full_list[,c("gene", "GO.Molecular.Function")]
  
  #get only non-empty ones (I just cannot analyze all the rest)
  if (gocat == "BP") BPGO <- BP[BP$GO.Biological.Process != "-",]
  if (gocat == "MF")BPGO <- BP[BP$GO.Molecular.Function != "-",]
  #get all GO terms in machine-readable way (in a list + without names, only numbers)
  
  if (gocat == "BP") GOs <- strsplit(BPGO$GO.Biological.Process, split = "| ", fixed = T)
  if (gocat == "MF") GOs <- strsplit(BPGO$GO.Molecular.Function, split = "| ", fixed = T)
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
  allRes <- GenTable(object = GOdata, f, topNodes = signif_at_0.001, numChar = 100)        
  write.csv(allRes, paste0("GO/", filename, "GO", upORdown, ".csv"))
  return(allRes)
}



# #Example code: get genes associated with a particular GO term
# allGO <- genesInTerm(GOdata)
# cdgenes <- allGO["GO:0046686"]
# allcdgenes <- full_list[cdgenes$`GO:0046686` %in% full_list$Sequence.Name,]
# decdgenes <- allcdgenes[allcdgenes$padj < 0.001 & allcdgenes$log2FoldChange > 3,]


#change here for another directory
#setwd("~/Documents/Paper1_stresses/data_tables/")
#FA_filename = "/homes/brauerei/polina/Documents/transcriptome_annotation/EveBCdTP_AnnotationTableFull.txt"
# setwd("~/Documents/Paper1_stresses/data_tables/by_eve")
# FA_filename = "../../../transcriptome_annotation/EveBCdTP1_cor/EveBCdTP1_cor_AnnotationTable.txt"
## For G. lacustris (assembly 2)
#setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/")
#FA_filename = "../../../transcriptome_annotation/GlaBCdTP1_cor2/GlaBCdTP1_cor2_AnnotationTable.txt"

setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/")
fa_filename = "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt"

setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced/")
fa_filename = "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt"

for (file in dir()) GOenrichment(filename = file, FA_filename = fa_filename)
for (file in dir()) GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")



setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")

for (file in dir()) {
  if(grepl("csv", file)) {
    fa_filename <- ifelse(startsWith(file, "Eve"), "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt",
                          ifelse(startsWith(file, "Ecy"), "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt",
                                 ifelse(startsWith(file, "Gla"), "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt", "")))
    GOenrichment(filename = file, FA_filename = fa_filename)
    GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")
  }
}


GOenrichment(filename = "./EvePhB03_annot.csv", FA_filename = "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt")
GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")

#setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/")

for (file in dir()) {
  fa_filename <- ifelse(startsWith(file, "Eve"), "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt",
                        ifelse(startsWith(file, "Ecy"), "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt",
                               ifelse(startsWith(file, "Gla"), "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt", "")))
  GOenrichment(filename = file, FA_filename = fa_filename)
  GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")
}
# 
# for (file in dir()) {
#   if(startsWith(file, "Gla")) {
#     fa_filename <- "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt"
#   GOenrichment(filename = file, FA_filename = fa_filename)
#   GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")}
# }



setwd("~/Documents/Paper1_stresses/data_tables/eve_deep/")
fa_filename = "~/Documents/transcriptome_annotation/EveBCdT_AnnotationTable.txt"

for (file in dir()) GOenrichment(filename = file, FA_filename = fa_filename)
for (file in dir()) GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")


setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/acetone_phenanthrene/")

for (file in dir()) {
  if(grepl("csv", file)) {
    fa_filename <- ifelse(startsWith(file, "Eve"), "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt",
                          ifelse(startsWith(file, "Ecy"), "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt",
                                 ifelse(startsWith(file, "Gla"), "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt", "")))
    GOenrichment(filename = file, FA_filename = fa_filename)
    GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")
  }
}

setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/acetone_phenanthrene/")

for (file in dir()) {
  if(grepl("csv", file)) {
    fa_filename <- ifelse(startsWith(file, "Eve"), "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt",
                          ifelse(startsWith(file, "Ecy"), "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt",
                                 ifelse(startsWith(file, "Gla"), "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt", "")))
    GOenrichment(filename = file, FA_filename = fa_filename)
    GOenrichment(filename = file, FA_filename = fa_filename, upORdown = "down")
  }
}

setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/test_sanity/")


####


















 
###########################################
  #rm(list = ls()) #ok, not now, but it's a good idea anyway

options(stringsAsFactors = F)
#options(device='x11')
#grDevices::X11.options(type='cairo')
options(bitmapType='cairo')

#    scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73",
#"#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
#black orange sky_blue green yellow blue vermillion reddish_purple

require(plyr)
require(reshape2)
require(ggplot2)
library(RColorBrewer)

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

######

#for now I need this one
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
                        paste0(substr(df$Term, 1, 20), "...", 
                               substr(df$Term, nchar(df$Term)-3, nchar(df$Term))),
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
  
  ggsave(paste0("../", condition, "_", upORdown, ".png"), width = width, height = height)
  ggsave(paste0("../", condition, "_", upORdown, ".svg"), width = width, height = height)
  
  print(p)
  return(p)
}

##manual switch here, which is bad. 
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO")
#setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced/GO")
upp <- gogo("LT10", "up")
downp <- gogo("LT10", "down")
upcd <- gogo("Cd", "up")
downcd <- gogo("Cd", "down")
upac <- gogo("PB", "up")
downac <- gogo("PB", "down")
upph <- gogo("Ph", "up")
downph <- gogo("Ph", "down")




setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/acetone_phenanthrene/GO")
upac <- gogo("PB", "up")
downac <- gogo("PB", "down")
upph <- gogo("Ph", "up")
downph <- gogo("Ph", "down")
upac <- gogo("PhB", "up")
downac <- gogo("PhB", "down")



setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/GO")
upp <- gogo("LT10", "up", suffix = "metazoa")
downp <- gogo("LT10", "down", suffix = "metazoa")
upcd <- gogo("Cd", "up", suffix = "metazoa")
downcd <- gogo("Cd", "down", suffix = "metazoa")
upac <- gogo("PB", "up", suffix = "metazoa")
downac <- gogo("PB", "down", suffix = "metazoa")
upph <- gogo("Ph", "up", suffix = "metazoa")
downph <- gogo("Ph", "down", suffix = "metazoa")


setwd("~/Documents/Paper1_stresses/data_tables/full_assemblies_each_by_own/acetone_phenanthrene/GO")
upac <- gogo("PB", "up", suffix = "metazoa")
downac <- gogo("PB", "down", suffix = "metazoa")
upph <- gogo("Ph", "up", suffix = "metazoa")
downph <- gogo("Ph", "down", suffix = "metazoa")
upphb <- gogo("PhB", "up", suffix = "metazoa")
downphb <- gogo("PhB", "down", suffix = "metazoa")



















####What if we take top500 instead of DE (accounting for effect size?)

GOenrichmentTop <- function(filename, sep = ",", upORdown = "up", gocat = "BP",
                         FA_filename = fa_filename) {
  #read full data set
  full_list <- read.csv(filename, sep = sep)
  ##rename the variable we'll need later for merging
  names(full_list)[3] <- "Sequence.Name" #2 for Gla, need to generate new ones
  #names(full_list)[7] <- "Best diamond hit"
  annotation <- read.delim(FA_filename)
  #merge tables
  withAnnotation <- merge(x = annotation, y=full_list[,c(2,3,4,7)], all.x=F, all.y=F , by="Sequence.Name")
  #get de genes
  #either upregs
  if (upORdown == "up") {
    de <- full_list[full_list$padj < 0.001,]
    de <- de[order(de$log2FoldChange, decreasing = T)[1:500], ]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "Sequence.Name"] 
  }
  #or downregs
  if (upORdown == "down") {
    de <- full_list[full_list$log2FoldChange < -3 & full_list$padj < 0.001,]
    de <- de[order(de$log2FoldChange)[1:500], ]
    #get rid of NAs & get names only
    de <- de[complete.cases(de$padj), "Sequence.Name"]
  }  
  
  
  #get GO terms
  if (gocat == "BP") BP <- withAnnotation[,c("Sequence.Name", "GO.Biological.Process")]
  if (gocat == "MF") BP <- withAnnotation[,c("Sequence.Name", "GO.Molecular.Function")]
  
  #get only non-empty ones (I just cannot analyze all the rest)
  if (gocat == "BP") BPGO <- BP[BP$GO.Biological.Process != "-",]
  if (gocat == "MF")BPGO <- BP[BP$GO.Molecular.Function != "-",]
  #get all GO terms in machine-readable way (in a list + without names, only numbers)
  
  if (gocat == "BP") GOs <- strsplit(BPGO$GO.Biological.Process, split = "| ", fixed = T)
  if (gocat == "MF") GOs <- strsplit(BPGO$GO.Molecular.Function, split = "| ", fixed = T)
  names(GOs) <- BPGO$Sequence.Name
  
  GOsTop <- lapply(X = GOs, function(x) gsub(" .*", "", x)) #remove human-readable name
  #get DE genes for the object
  DElist <- factor(as.integer(names(GOsTop) %in% de))
  
  #exit the function if the DE list is empty
  if (length(levels(DElist)) <= 1) return("No DE genes to look at")
  
  
  names(DElist) <- names(GOsTop)
  #construct a GOdat object (topGO package)
  GOdata <- new("topGOdata", ontology = gocat, allGenes = DElist, annot = annFUN.gene2GO, gene2GO = GOsTop)
  f <- runTest(GOdata, algorithm = "elim", statistic = "fisher") 
  #from the manual: We like to interpret the p-values returned by these methods as corrected or not affected by multiple testing.    
  signif_at_0.001 <- sum(f@score < 0.001)
  allRes <- GenTable(object = GOdata, f, topNodes = signif_at_0.001, numChar = 100)        
  write.csv(allRes, paste0("GOtop500/", filename, "GO", upORdown, ".csv"))
  return(allRes)
}


setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")

for (file in dir()) {
  if(grepl("csv", file)) {
    fa_filename <- ifelse(startsWith(file, "Eve"), "~/Documents/transcriptome_annotation/EveBCdTP1_cor_AnnotationTable.txt",
                          ifelse(startsWith(file, "Ecy"), "~/Documents/transcriptome_annotation/EcyBCdTP1_cor_AnnotationTable.txt",
                                 ifelse(startsWith(file, "Gla"), "~/Documents/transcriptome_annotation/GlaBCdTP1_cor2_AnnotationTable.txt", "")))
    GOenrichmentTop(filename = file, FA_filename = fa_filename)
    GOenrichmentTop(filename = file, FA_filename = fa_filename, upORdown = "down")
  }

}


##manual switch here, which is bad. 
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO")
#setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced/GO")
upp <- gogo("LT10", "up")
downp <- gogo("LT10", "down")
upcd <- gogo("Cd", "up")
downcd <- gogo("Cd", "down")
upac <- gogo("PB", "up")
downac <- gogo("PB", "down")
upph <- gogo("Ph", "up")
downph <- gogo("Ph", "down")
upph <- gogo("PhB", "up")
downph <- gogo("PhB", "down")

setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO")
#setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced/GO")
upp <- gogo("LT10", "up")
downp <- gogo("LT10", "down")
upcd <- gogo("Cd", "up")
downcd <- gogo("Cd", "down")
upac <- gogo("PB", "up")
downac <- gogo("PB", "down")
upph <- gogo("Ph", "up")
downph <- gogo("Ph", "down")



























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

setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/GO")


gogots(species = "Eve")
gogots(species = "Ecy")
gogots(species = "Gla")

gogots(species = "Eve", upORdown = "down")
gogots(species = "Ecy", upORdown = "down")
gogots(species = "Gla", upORdown = "down")
