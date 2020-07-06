options(stringsAsFactors = F)
library(UpSetR)
library(ggplot2)

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


    ## read files
    ## salmon / DESeq2
    sd <- read.delim(paths[1])

    #select DE genes...
    sd_DE <- row.names(selectDEd(sd))

    ## intersections at least
    ## upSet figures
    lI <- list(salmon.DESeq2 = sd_DE, salmon.edgeR = se_DE, 
               rsem.DESeq2 = rd_DE, rsem.edgeR = re_DE)
    upset(fromList(lI), nintersects = NA)
    ## get the list of common elements
    common_list <- Reduce(intersect, lI)
    
    sd_common <- sd[row.names(sd) %in% common_list, ]
    sd_common <- sd_common[order(row.names(sd_common)), ]
  
    ## sanity check!
    print(row.names(sd_common) == row.names(re_common))
    
    result <- sd.common
    
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



# #install.packages("openxlsx")
# library(openxlsx)
# write.xlsx(x = list("Ecy_ace_3h" = ecy_a03, "Ecy_ace_24h" = ecy_a24, "Ecy_phe_3h" = ecy_p03, "Ecy_phe_24h"= ecy_p24, 
#                     "Ecy_ace+phe_3h" = ecy_ap03, "Ecy_ace+phe_24h" = ecy_ap24, 
#                     "Eve_ace_24h" = eve_a24, "Eve_phe_3h" = eve_p03, "Eve_phe_24h"= eve_p24, 
#                     "Eve_ace+phe_3h" = eve_ap03, "Eve_ace+phe_24h" = eve_ap24, 
#                     "Gla_ace_3h" = gla_a03, "Gla_ace_24h" = gla_p24, "Gla_phe_3h" = gla_p03, "Gla_phe_24h"= gla_p24, 
#                     "Gla_ace+phe_24h" = gla_ap24),
#            file = "DE.xlsx")


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