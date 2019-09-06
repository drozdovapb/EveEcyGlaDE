library(ggplot2)
options(stringsAsFactors = F)
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#correlation?

plotcor <- function(thisspecies = "Eve") {
  EveT <- read.csv(paste0("./", thisspecies, "LT1024_annot.csv"))
  EveCd <- read.csv(paste0("./", thisspecies, "Cd24_annot.csv"))
  
  EveTDE <- read.csv(paste0("./DE/", thisspecies, "LT1024_annot_de.csv"))
  EveCdDE <- read.csv(paste0("./DE/", thisspecies, "Cd24_annot_de.csv"))
  
  
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thiscolor <- colors[which(species == thisspecies)]
  
  ## and now select only DE genes!!
  
  degenes <- unique(c(EveTDE$gene, EveCdDE$gene))
  
  merged <- merge(EveT, EveCd, by = "gene")
  #geom_point(data = tbl[tbl$ishsp70,], shape = 8, col = "#cc79a7", alpha = 1) #pink
  
  mergedDE <- merged[merged$gene %in% degenes, ]
  write.csv(mergedDE, paste0(thisspecies, "_T_vs_Cd.csv"))
  
  
  xmin <- min(mergedDE$log2FoldChange.x)
  xmax <- max(mergedDE$log2FoldChange.x)
  ymin <- min(mergedDE$log2FoldChange.y)
  ymax <- max(mergedDE$log2FoldChange.y)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
  p <- ggplot(data = mergedDE, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(alpha = .1, col = "black") + geom_smooth(method = "lm", col = "darkgray") + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Elevated temperature, 24 hours") + ylab("Cadmium, 24 hours") + 
    theme_bw(base_size = 16) + xlim(fmin, fmax) + ylim(fmin, fmax)
  
  
  
  p0 <- ggplot(data = mergedDE, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    #geom_point(alpha = .1, col = "black") + geom_smooth(method = "lm", col = "darkgray") + 
    geom_text(label = mergedDE$best.nr.hit.diamond.x, size = 2) + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Elevated temperature, 24 hours") + ylab("Cadmium, 24 hours") + 
    theme_bw(base_size = 16) + xlim(xmin, xmax) + ylim(ymin, ymax)
  
  #and now get 
  is.translation <- grepl("translation", mergedDE$GO.Biological.Process.x)
    #grepl("GO:0006414 translational elongation", mergedDE$GO.Biological.Process.x)
  is.chitin <- grepl("chitin metabolic process", mergedDE$GO.Biological.Process.x)
  #is.metabo <- grepl("carbohydrate metabolic process", mergedDE$GO.Biological.Process.x)
  is.stress <- grepl("response to stress", mergedDE$GO.Biological.Process.x) | grepl("response to heat", mergedDE$GO.Biological.Process.x)
  is.proteo <- grepl("proteolysis", mergedDE$GO.Biological.Process.x)
  is.ubiqui <- grepl("ubiquitin", mergedDE$GO.Biological.Process.x)
  
  p1 <- p + geom_point(data = mergedDE[is.translation,], shape = 5, color = thiscolor) + 
    geom_point(data = mergedDE[is.chitin, ], shape = 6, color = thiscolor) + 
    geom_point(data = mergedDE[is.stress, ], shape = 8, color = thiscolor) + 
    geom_point(data = mergedDE[is.proteo, ], shape = 10, color = thiscolor)

  p2 <- p + geom_text(data = mergedDE[is.translation,], label = "T", color = thiscolor) + 
    geom_text(data = mergedDE[is.chitin, ], label = "C", color = thiscolor) + 
    geom_text(data = mergedDE[is.stress, ], label = "S", color = thiscolor) + 
    geom_text(data = mergedDE[is.proteo, ], label = "P", color = thiscolor) + 
    geom_text(data = mergedDE[is.ubiqui, ], label = "U", color = thiscolor)
  
    
  print(p2)
  #print(p)
  
  cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")
  summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))
  
  return(p2)
}

pv <- plotcor()
pc <- plotcor(thisspecies = "Ecy")
pl <- plotcor(thisspecies = "Gla")


ggsave("Eve_corr.svg", pv, width = 8, height = 8)
ggsave("Ecy_corr.svg", pc, width = 4, height = 4)
ggsave("Gla_corr.svg", pl, width = 4, height = 4)
