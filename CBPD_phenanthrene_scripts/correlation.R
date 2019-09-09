library(ggplot2)
options(stringsAsFactors = F)
setwd("/media/drozdovapb/big/Research/Projects/DE/texts/Paper1_stresses/data_tables/reduced_each_by_own/")
#correlation?
defold <- "/media/drozdovapb/big/Research/Projects/DE/texts/Paper1_stresses/acetone_phenanthrene_story/csv/"


plotcor <- function(thisspecies = "Eve") {
  EveT <- read.csv(paste0("./", thisspecies, "PB24_annot.csv"))
  EveCd <- read.csv(paste0("./", thisspecies, "Ph24_annot.csv"))
  
  EveTDE <- read.csv(paste0(defold, thisspecies, "PB24_annot_de.csv"))
  EveCdDE <- read.csv(paste0(defold, thisspecies, "Ph24_annot_de.csv"))
  
  
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thiscolor <- colors[which(species == thisspecies)]
  
  ## and now select only DE genes!!
  
  degenes <- unique(c(EveTDE$gene, EveCdDE$gene))
  
  merged <- merge(EveT, EveCd, by = "gene")
  #geom_point(data = tbl[tbl$ishsp70,], shape = 8, col = "#cc79a7", alpha = 1) #pink
  
  mergedDE <- merged[merged$gene %in% degenes, ]
  #write.csv(mergedDE, paste0(thisspecies, "_T_vs_Cd.csv"))
  
  
  xmin <- min(mergedDE$log2FoldChange.x)
  xmax <- max(mergedDE$log2FoldChange.x)
  ymin <- min(mergedDE$log2FoldChange.y)
  ymax <- max(mergedDE$log2FoldChange.y)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
  fit <- lm(merged$log2FoldChange.y~merged$log2FoldChange.x)
  
  p <- ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(size = 0.5, alpha = .2, col = "black") + geom_smooth(method = "lm", col = "darkgray") + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Acetone, 24 hours") + ylab("Phenanthrene, 24 hours") + 
    theme_bw(base_size = 12) + xlim(fmin, fmax) + ylim(fmin, fmax) + 
    labs(aes(size = 4), title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + 
    coord_fixed()

  
  p
  
  print(cor.test(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson"))
  print(cor.test(merged$log2FoldChange.y, merged$log2FoldChange.x, method = "pearson"))
  #print(summary(lm(merged$log2FoldChange.y ~ merged$log2FoldChange.x)))
  
  print(p)  
  return(p)
}

pv <- plotcor() 
savedir <- "/media/drozdovapb/big/Research/Projects/DE/texts/Paper1_stresses/acetone_phenanthrene_story/figures/correlation/"
ggsave(paste(savedir, "verr_all.svg"), width = 10)
ggsave(paste0(savedir, "Eve_all.png"), width = 10)

pc <- plotcor(thisspecies = "Ecy")
ggsave(paste0(savedir, "Ecy_all.png"), width = 10)
pl <- plotcor(thisspecies = "Gla")
ggsave(paste0(savedir, "Gla_all.png"), width = 10)


plotcorDE <- function(thisspecies = "Eve") {
  EveT <- read.csv(paste0("./", thisspecies, "PB24_annot.csv"))
  EveCd <- read.csv(paste0("./", thisspecies, "Ph24_annot.csv"))
  
  EveTDE <- read.csv(paste0(defold, thisspecies, "PB24_annot_de.csv"))
  EveCdDE <- read.csv(paste0(defold, thisspecies, "Ph24_annot_de.csv"))
  
  
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thiscolor <- colors[which(species == thisspecies)]
  
  ## and now select only DE genes!!
  
  degenes <- unique(c(EveTDE$gene, EveCdDE$gene))
  
  merged <- merge(EveT, EveCd, by = "gene")
  #geom_point(data = tbl[tbl$ishsp70,], shape = 8, col = "#cc79a7", alpha = 1) #pink
  
  mergedDE <- merged[merged$gene %in% degenes, ]
  #write.csv(mergedDE, paste0(thisspecies, "_T_vs_Cd.csv"))
  
  
  xmin <- min(mergedDE$log2FoldChange.x)
  xmax <- max(mergedDE$log2FoldChange.x)
  ymin <- min(mergedDE$log2FoldChange.y)
  ymax <- max(mergedDE$log2FoldChange.y)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
  fit <- lm(mergedDE$log2FoldChange.y~mergedDE$log2FoldChange.x)
  print(summary(fit))
  
  p <- ggplot(data = mergedDE, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(size = .5, alpha = .5, col = "black") + geom_smooth(method = "lm", col = "darkgray") + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Acetone, 24 hours") + ylab("Phenanthrene, 24 hours") + 
    theme_bw(base_size = 12) + xlim(fmin, fmax) + ylim(fmin, fmax) + 
    #labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                   "Intercept =",signif(fit$coef[[1]],5 ),
    #                   " Slope =",signif(fit$coef[[2]], 5),
    #                   " P =",signif(summary(fit)$coef[2,4], 5))) + 
    coord_fixed()
  
  #print(summary(lm(mergedDE$log2FoldChange.x ~ mergedDE$log2FoldChange.y)))
  print(cor.test(mergedDE$log2FoldChange.x, mergedDE$log2FoldChange.y, method = "pearson"))
  print(cor.test(mergedDE$log2FoldChange.y, mergedDE$log2FoldChange.x, method = "pearson"))
  
  
  print(p)  
  return(p)
}

pv <- plotcorDE()
ggsave(paste0(savedir,"EveDE.svg"), pv, width = 10)
pc <- plotcorDE(thisspecies = "Ecy")
ggsave(paste0(savedir,"EcyDE.svg"), pc, width = 10)
pv <- plotcorDE(thisspecies = "Gla")
ggsave(paste0(savedir,"GlaDE.svg"), pl, width = 10)
ggsave(paste0(savedir,"GlaDE.png"), pl, width = 10)





plotcorDE2 <- function(thisspecies = "Eve", thiscolor = "red4") {
  EveT <- read.csv(paste0("./", thisspecies, "PB24_annot.csv"))
  EveCd <- read.csv(paste0("./", thisspecies, "Ph24_annot.csv"))
  
  EveTDE <- read.csv(paste0(defold, thisspecies, "PB24_annot_de.csv"))
  EveCdDE <- read.csv(paste0(defold, thisspecies, "Ph24_annot_de.csv"))
  
  
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thiscolor <- colors[which(species == thisspecies)]
  
  ## and now select only DE genes!!
  
  degenes <- unique(c(EveTDE$gene, EveCdDE$gene))
  
  merged <- merge(EveT, EveCd, by = "gene")
  #geom_point(data = tbl[tbl$ishsp70,], shape = 8, col = "#cc79a7", alpha = 1) #pink
  
  mergedDE <- merged[merged$gene %in% degenes, ]
  #write.csv(mergedDE, paste0(thisspecies, "_T_vs_Cd.csv"))
  
  
  xmin <- min(mergedDE$log2FoldChange.x)
  xmax <- max(mergedDE$log2FoldChange.x)
  ymin <- min(mergedDE$log2FoldChange.y)
  ymax <- max(mergedDE$log2FoldChange.y)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
  fit <- lm(mergedDE$log2FoldChange.y~mergedDE$log2FoldChange.x)
  print(summary(fit))
  
  p <- ggplot(data = mergedDE, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(size = .5, alpha = .5, col = thiscolor) + geom_smooth(method = "lm", col = thiscolor) + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Acetone, 24 hours") + ylab("Phenanthrene, 24 hours") + 
    theme_bw(base_size = 12) + xlim(fmin, fmax) + ylim(fmin, fmax) + 
    #labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                   "Intercept =",signif(fit$coef[[1]],5 ),
    #                   " Slope =",signif(fit$coef[[2]], 5),
    #                   " P =",signif(summary(fit)$coef[2,4], 5))) + 
    coord_fixed()
  
  #print(summary(lm(mergedDE$log2FoldChange.x ~ mergedDE$log2FoldChange.y)))
  print(cor.test(mergedDE$log2FoldChange.x, mergedDE$log2FoldChange.y, method = "pearson"))
  print(cor.test(mergedDE$log2FoldChange.y, mergedDE$log2FoldChange.x, method = "pearson"))
  
  
  print(p)  
  return(p)
}


savedir <- "/media/drozdovapb/big/Research/Projects/DE/texts/Paper1_stresses/acetone_phenanthrene_story/CBPD_submit/temp/"

pv <- plotcorDE2(thisspecies = "Eve", thiscolor = "#4087AF") 
ggsave(paste0(savedir,"EveDE_abstr.svg"), pv, width = 2.5, height = 2.5)
pc <- plotcorDE2(thisspecies = "Ecy", thiscolor = "#007656")
ggsave(paste0(savedir,"EcyDE_abstr.svg"), pc, width = 2.5, height = 2.5)
pl <- plotcorDE2(thisspecies = "Gla", thiscolor = "#D55E00")
ggsave(paste0(savedir,"GlaDE_abstr.svg"), pl, width = 2.5, height = 2.5)


#"#4087AF", "#007656", "#D55E00"