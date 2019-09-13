options(stringsAsFactors = F)
library(UpSetR)
library(ggplot2)

## re = rsem /edger;
## rd = rsem / deseq2; 
## se = salmon / edgr;
## sd = salmon / deseq2

## later make it a function
comparison <- "PB24_vs_B24h"
species <- "Ecy"

## a small function to extract DE genes
threshold_p <- 0.05
selectDEd <- function(tbl) {
  tbl[(tbl$log2FoldChange > 1 | tbl$log2FoldChange < -1) & tbl$padj < threshold_p, ]
}
selectDEe <- function(tbl) {
  tbl[(tbl$logFC > 1 | tbl$logFC < -1) & tbl$FDR < threshold_p, ]
}


get_common_list <- function(comparison, species) {
      ## I have 4 files to compare, just need to read them all and combine
    est_method <- c("salmon", "salmon", "bowtie_rsem", "bowtie_rsem")
    de_method <- c("DESeq2", "edgeR", "DESeq2", "edgeR")
    ## paths to 4 files
    paths <- paste0(est_method, "/", species, ".isoform.counts.matrix.", 
                    comparison, ".", de_method, ".DE_results")
    ## read files
    ## salmon / DESeq2
    sd <- read.delim(paths[1])
    ## salmon / edgeR
    se <- read.delim(paths[2])
    ## bowtie_rsem / DESeq2
    rd <- read.delim(paths[3])
    ## bowtie_rsem / edgeR
    re <- read.delim(paths[4])
  
    #select DE genes...
    sd_DE <- row.names(selectDEd(sd))
    se_DE <- row.names(selectDEe(se))
    rd_DE <- row.names(selectDEd(rd))
    re_DE <- row.names(selectDEe(re))
  
    ## intersections at least
    ## upSet figures
    lI <- list(salmon.DESeq2 = sd_DE, salmon.edgeR = se_DE, 
               rsem.DESeq2 = rd_DE, rsem.edgeR = re_DE)
    upset(fromList(lI), nintersects = NA)
    ## get the list of common elements
    common_list <- Reduce(intersect, lI)
    
    sd_common <- sd[row.names(sd) %in% common_list, ]
    sd_common <- sd_common[order(row.names(sd_common)), ]
    se_common <- se[row.names(se) %in% common_list, ]
    se_common <- se_common[order(row.names(se_common)), ]
    rd_common <- rd[row.names(rd) %in% common_list, ]
    rd_common <- rd_common[order(row.names(rd_common)), ]
    re_common <- re[row.names(re) %in% common_list, ]
    re_common <- re_common[order(row.names(re_common)), ]
    
    ## sanity check!
    print(row.names(sd_common) == row.names(re_common))
    
    result <- data.frame(Gene = row.names(sd_common),
                         salmon.DESeq2.FC = sd_common$log2FoldChange, 
                         salmon.DESeq2.padj = sd_common$padj,
                         salmon.edgeR.FC = se_common$logFC,
                         salmon.edgeR.FDR = se_common$FDR,
                         rsem.DESeq2.FC = rd_common$log2FoldChange, 
                         rsem.DESeq2.padj = rd_common$padj,
                         rsem.edgeR.FC = re_common$logFC, 
                         rsem.edgeR.FDR = re_common$FDR)
    
## if it works, then add p-values and so on!
    # the line below: write just once
    diamond <- read.delim(paste0("~/Research/Projects/DE/annotation/",
                                 species, "BCdTP1_cor.diamond.tsv"), head = F)
    ## first check if the '%in%' method works
    ##result$best.match.diamond <- diamond[diamond$V1 %in% common_list, "V1"]
    ##sort(result$best.match.diamond) == sort(result$Gene)
        ## nope, it doesn't. The order is wrong.
    ##result$check.diamond <- diamond[diamond$V1 %in% common_list, "V1"]    
    ##result$best.match.diamond <- sapply(common_list, 
      ##                          function(x) diamond[diamond$V1 %in% x, "V1"])
    ##sort(result$best.match.diamond) == sort(result$Gene)
    ## now it works!
    result$best.match.diamond <- sapply(as.character(result$Gene), 
                              function(x) diamond[diamond$V1 %in% x, "V14"])
    
    fa <- read.delim(paste0("~/Research/Projects/DE/annotation/",
                                 species, "BCdTP1_cor_AnnotationTable.txt"), head = T)
    result$best.match.to.nr <- sapply(as.character(result$Gene), 
                              function(x) fa[fa$Sequence.Name %in% x, "best.hit.to.nr"])
    result$GO.Biological.Process <- sapply(as.character(result$Gene), 
                      function(x) fa[fa$Sequence.Name %in% x, "GO.Biological.Process"])
    result$GO.Molecular.Function <- sapply(as.character(result$Gene), 
                     function(x) fa[fa$Sequence.Name %in% x, "GO.Molecular.Function"])
    result$GO.Cellular.Component <- sapply(as.character(result$Gene), 
                             function(x) fa[fa$Sequence.Name %in% x, "GO.Cellular.Component"])
    
    write.csv(result, paste0(comparison, ".", species, "common_DE.csv"))
    return(result)
}

## Get all common lists 
ecy_a03 <- get_common_list(comparison = "PB03_vs_B03", species = "Ecy")
ecy_a24 <- get_common_list("PB24_vs_B24h", "Ecy")
ecy_p03 <- get_common_list("Ph03_vs_PB03", "Ecy")
ecy_p24 <- get_common_list("Ph24_vs_PB24", "Ecy")
ecy_ap03 <- get_common_list("Ph03_vs_B03", "Ecy")
ecy_ap24 <- get_common_list("Ph24_vs_B24h", "Ecy")

eve_a03 <- get_common_list(comparison = "PB03_vs_B03", species = "Eve")
eve_a24 <- get_common_list("PB24_vs_B24h", "Eve")
eve_p03 <- get_common_list("Ph03_vs_PB03", "Eve")
eve_p24 <- get_common_list("Ph24_vs_PB24", "Eve")
eve_ap03 <- get_common_list("Ph03_vs_B03", "Eve")
eve_ap24 <- get_common_list("Ph24_vs_B24h", "Eve")

gla_a03 <- get_common_list(comparison = "PB03_vs_B03", species = "Gla")
gla_a24 <- get_common_list("PB24_vs_B24h", "Gla")
gla_p03 <- get_common_list("Ph03_vs_PB03", "Gla")
gla_p24 <- get_common_list("Ph24_vs_PB24", "Gla")
gla_ap03 <- get_common_list("Ph03_vs_B03", "Gla")
gla_ap24 <- get_common_list("Ph24_vs_B24h", "Gla")

## should be a possible option but doesn't work here due to some java issues
##library(xlsx)
#index = 1
#for (file in dir()) {
#  if (grepl(".csv", file)) {
#    temp <- read.csv(file)
#    if (nrow(temp) > 0) {
#      ifappend = index > 1
#      write.xlsx(temp, "DE.xlsx", append = ifappend, sheetName = file)
#      index <- index + 1
#    }
#  }
#}

#install.packages("openxlsx")
library(openxlsx)
write.xlsx(x = list("Ecy_ace_3h" = ecy_a03, "Ecy_ace_24h" = ecy_a24, "Ecy_phe_3h" = ecy_p03, "Ecy_phe_24h"= ecy_p24, 
                    "Ecy_ace+phe_3h" = ecy_ap03, "Ecy_ace+phe_24h" = ecy_ap24, 
                    "Eve_ace_24h" = eve_a24, "Eve_phe_3h" = eve_p03, "Eve_phe_24h"= eve_p24, 
                    "Eve_ace+phe_3h" = eve_ap03, "Eve_ace+phe_24h" = eve_ap24, 
                    "Gla_ace_3h" = gla_a03, "Gla_ace_24h" = gla_p24, "Gla_phe_3h" = gla_p03, "Gla_phe_24h"= gla_p24, 
                    "Gla_ace+phe_24h" = gla_ap24),
           file = "DE.xlsx")

## I'm afraid but need to do it. What's the correlation?


correlation <- function(ace, phe, aceDE, pheDE) {
  
  ace$Gene <- row.names(ace)
  phe$Gene <- row.names(phe)
  merged <- merge(ace, phe, by = "Gene", all = T)
  
  DElist <- unique(c(aceDE$Gene, pheDE$Gene))
  mergedDE <- merged[merged$Gene %in% DElist, ]
  
  xmin <- min(mergedDE$log2FoldChange.x, na.rm = T)
  xmax <- max(mergedDE$log2FoldChange.x,  na.rm = T)
  ymin <- min(mergedDE$log2FoldChange.y,  na.rm = T)
  ymax <- max(mergedDE$log2FoldChange.y,  na.rm = T)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
  fit <- lm(mergedDE$log2FoldChange.y~mergedDE$log2FoldChange.x)
  print(summary(fit))
  
  p <- ggplot(data = mergedDE, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(size = .5, alpha = .5, col = "black") + geom_smooth(method = "lm", col = "gray") + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Solvent control, 24 hours") + ylab("Phenanthrene, 24 hours") + 
    theme_bw(base_size = 12) + xlim(fmin, fmax) + ylim(fmin, fmax) + 
    #labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                   "Intercept =",signif(fit$coef[[1]],5 ),
    #                   " Slope =",signif(fit$coef[[2]], 5),
    #                   " P =",signif(summary(fit)$coef[2,4], 5))) + 
    coord_fixed()
  
  #print(summary(lm(mergedDE$log2FoldChange.x ~ mergedDE$log2FoldChange.y)))
  print(cor.test(mergedDE$log2FoldChange.x, mergedDE$log2FoldChange.y, method = "pearson"))
  return (p)
  
}


species = "Eve"
eve_ace <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                   "PB24_vs_B24h", ".", "DESeq2", ".DE_results"))
eve_phe <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                  "Ph24_vs_PB24", ".", "DESeq2", ".DE_results"))
correlation(ace = eve_ace, phe = eve_phe, aceDE = eve_a24, pheDE = eve_p24)
ggsave("Eve_corr.svg", width = 3, height = 3)


species <- "Ecy"
ecy_ace <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                             "PB24_vs_B24h", ".", "DESeq2", ".DE_results"))
ecy_phe <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                             "Ph24_vs_PB24", ".", "DESeq2", ".DE_results"))
correlation(ace = ecy_ace, phe = ecy_phe, aceDE = ecy_a24, pheDE = ecy_p24)
ggsave("Ecy_corr.svg", width = 3, height = 3)


species <- "Gla"
gla_ace <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                             "PB24_vs_B24h", ".", "DESeq2", ".DE_results"))
gla_phe <- read.delim(paste0("salmon/", species, ".isoform.counts.matrix.", 
                             "Ph24_vs_PB24", ".", "DESeq2", ".DE_results"))
correlation(ace = gla_ace, phe = gla_phe, aceDE = gla_a24, pheDE = gla_p24)
ggsave("Gla_corr.svg", width = 4, height = 4)
ggsave("Gla_corr.png", width = 3, height = 3)



## figure 1 (numer of DE genes)
## eg 
## table(ecy_ap24$salmon.edgeR.FC > 0)


numberDE <- read.csv("./ace_phe_DE.csv")
numberDE$samplegroup <- paste(numberDE$Species, numberDE$Time)


ace <- numberDE[numberDE$Treatment=="ace", ]
phe <- numberDE[numberDE$Treatment=="phe", ]
ace_phe <- numberDE[numberDE$Treatment=="ace+phe", ]


ggplot(ace, aes(x = X, fill=Species)) + 
  geom_hline(yintercept = c(-80, -60, -40, -20, 20, 40, 60, 80), col = "grey") +
  geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
  geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
  ylab("Number of DE transcripts") + 
  geom_hline(yintercept = 0, col = "grey") + 
  ylim(-100, 100) 

ggsave("ace.svg",width=10, heigh=6, units = "cm")  

ggplot(phe, aes(x = X, fill=Species)) + 
  geom_hline(yintercept = c(-140, -100, -60, -20, 20, 60, 100, 140, 180), col = "grey") +
  geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
  geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
 ylab("Number of DE transcripts") + 
  geom_hline(yintercept = 0, col = "grey") + 
  ylim(-150, 200) 

ggsave("phe.svg",width=10, heigh=6, units = "cm")  


ggplot(ace_phe, aes(x = X, fill=Species)) + 
  geom_hline(yintercept = c(-60, -40, -20, 20, 40, 60), col = "grey") +
  geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
  geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
  ylab("Number of DE transcripts") + 
  geom_hline(yintercept = 0, col = "grey") + 
  ylim(-75, 75) 

ggsave("ace_phe.svg",width=10, heigh=6, units = "cm")  
