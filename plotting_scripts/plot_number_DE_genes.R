#############################################################3
###environment
options(stringsAsFactors = F)


###The main function
###Plot 6 pairs of bars, i.e. 3 species x 2 time points

getnumberDE <- function(conditions) { #conditions is a vector of length 2, e.g. c("LT1003", "LT1024")
  #create an empty data frame of the necessary dimensions
  numberDE <-data.frame(condition = rep(conditions, 3), 
                         species = c(rep("Eve", 2), rep("Ecy", 2), rep("Gla", 2))) 
  ##Some hard-coded examples
  # numberDE <-data.frame(condition = rep(c("LT1003", "LT1024"), 3), 
  #                       species = c(rep("Eve", 2), rep("Ecy", 2), rep("Gla", 2))) 
  # numberDE <-data.frame(condition = rep(c("Cd03", "Cd24"), 3), 
  #                       species = c(rep("Eve", 2), rep("Ecy", 2), rep("Gla", 2))) 

  ###remove replicates @was necessart in some case... 
  #numberDE <- numberDE[!duplicated(numberDE), ]
  #create columns
  numberDE$upAll <- NA
  numberDE$upAllMost <- NA
  numberDE$downAll <- NA
  numberDE$downAllMost <- NA
  #factors for plotting
  numberDE$samplegroup <- paste0(numberDE$species, numberDE$condition)

  ##populate the data frame (change the directory if needed)
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    species <- thisrow[2]
    ### in case you want everything by one reference, please use the next line:
    #(wdir <- paste0("/run/media/polina/Elements/transcriptome/DE/9-reduced/", tolower(species)))
    wdir <- "~/Documents/Paper1_stresses/data_tables/reduced_each_by_own" 
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot.csv"))
    numberDE$All[i] <- nrow(temp)
    numberDE$upAll[i] <- sum(temp$log2FoldChange > 1 & temp$padj < 0.001)
    numberDE$upAllMost[i] <- sum(temp$log2FoldChange > 3 & temp$padj < 0.001)
    numberDE$downAll[i] <- -sum(temp$log2FoldChange < -1 & temp$padj < 0.001)
    numberDE$downAllMost[i] <- -sum(temp$log2FoldChange < -3 & temp$padj < 0.001)
  }
  ##reverse factor levels for correct down=>up plotting
  numberDE$samplegroup <- factor(numberDE$samplegroup, levels = rev(numberDE$samplegroup)) #mainly for plotting
  pa <- ggplot(numberDE, aes(x = samplegroup, fill=species, col = species)) + 
    geom_bar(aes(y=upAll), stat="identity", position="dodge", alpha = 0.4) + 
    geom_bar(aes(y=upAllMost, col = species), stat="identity", position="dodge", alpha = 1) + 
    geom_bar(aes(y=downAll), stat="identity", position="dodge", alpha = 0.4) +
    geom_bar(aes(y=downAllMost, col = species), stat="identity", position="dodge", alpha = 1) + 
    theme_classic(base_size = 16) + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) +
    scale_color_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    ylab("number of DE genes, two-fold change") + 
    geom_hline(yintercept = 0, col = "white") #+ #ylim(-100, 100) +
  
  print(pa)
  
  ggsave(paste0("~/Documents/Paper1_stresses/multipanel_figures/", conditions[1], conditions[2], "_nDE.svg"), 
         pa, width = 150, units = "mm")
  
}

###############################################################################
###Main part
###Example usage
getnumberDE(c("LT1003", "LT1024"))
