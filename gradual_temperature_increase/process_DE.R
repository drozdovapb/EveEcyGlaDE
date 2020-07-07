options(stringsAsFactors = F)
#library(UpSetR)
library(ggplot2)
library(reshape2)

##setwd("~/Downloads/Primers/")


## later make it a function
#comparison <- "PB24_vs_B24h"
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
Ecy12 <- read.delim("Ecy.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results")
  Ecy12$Gene <- row.names(Ecy12)
Ecy18 <- read.delim("Ecy.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results")
  Ecy18$Gene <- row.names(Ecy18)
Ecy24 <- read.delim("Ecy.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results")
  Ecy24$Gene <- row.names(Ecy24)

Eve12 <- read.delim("Eve.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results")
  Eve12$Gene <- row.names(Eve12)
Eve18 <- read.delim("Eve.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results")
  Eve18$Gene <- row.names(Eve18)
Eve24 <- read.delim("Eve.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results")
  Eve24$Gene <- row.names(Eve24)

  
Gla12 <- read.delim("Gla.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results")
  Gla12$Gene <- row.names(Gla12)
Gla18 <- read.delim("Gla.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results")
  Gla18$Gene <- row.names(Gla18)
Gla24 <- read.delim("Gla.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results")
  Gla24$Gene <- row.names(Gla24)

  
    # #select DE genes...
    # sd_DE <- row.names(selectDEd(Ecy12))
  ### make full DE analysis

      # result <- Ecy24
      # result$Gene <- row.names(result)
      # 
    diamondEcy <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                 "Ecy", "BCdTP1_ani.diamond.tsv"), head = F)
    diamondEve <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                    "Eve", "BCdTP1_ani.diamond.tsv"), head = F)
    diamondGla <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                    "Gla", "BCdTP1_ani.diamond.tsv"), head = F)
    
    RPsEcy <- diamondEcy[grep("ribosomal pr", diamondEcy$V14), ]
    RPsEve <- diamondEve[grep("ribosomal pr", diamondEve$V14), ]
    RPsGla <- diamondGla[grep("ribosomal pr", diamondGla$V14), ]
    ## add the option with the GO term
    
    faEcy <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                                             "Ecy", "BCdTP1_cor_AnnotationTable.txt"), head = T)
    grep("GO:0006412", faEcy$GO.Biological.Process)
    
    
    EcyT12.RPs <- Ecy12[Ecy12$Gene %in% RPsEcy$V1, ]
    EcyT18.RPs <- Ecy18[Ecy18$Gene %in% RPsEcy$V1, ]
    EcyT24.RPs <- Ecy24[Ecy24$Gene %in% RPsEcy$V1, ]
    #only common genes
    commonEcy.RPs <- Reduce(x=list(EcyT12.RPs$Gene, EcyT18.RPs$Gene, EcyT24.RPs$Gene), f="intersect")
    EcyT12.RPs <- EcyT12.RPs[EcyT12.RPs$Gene %in% commonEcy.RPs, ]
    EcyT18.RPs <- EcyT18.RPs[EcyT18.RPs$Gene %in% commonEcy.RPs, ]
    EcyT24.RPs <- EcyT24.RPs[EcyT24.RPs$Gene %in% commonEcy.RPs, ]
    
    
    
    EveT12.RPs <- Eve12[Eve12$Gene %in% RPsEve$V1, ]
    EveT18.RPs <- Eve18[Eve18$Gene %in% RPsEve$V1, ]
    EveT24.RPs <- Eve24[Eve24$Gene %in% RPsEve$V1, ]
    #only common genes
    commonEve.RPs <- Reduce(x=list(EveT12.RPs$Gene, EveT18.RPs$Gene, EveT24.RPs$Gene), f="intersect")
    EveT12.RPs <- EveT12.RPs[EveT12.RPs$Gene %in% commonEve.RPs, ]
    EveT18.RPs <- EveT18.RPs[EveT18.RPs$Gene %in% commonEve.RPs, ]
    EveT24.RPs <- EveT24.RPs[EveT24.RPs$Gene %in% commonEve.RPs, ]
    
        
    GlaT12.RPs <- Gla12[Gla12$Gene %in% RPsGla$V1, ]
    GlaT18.RPs <- Gla18[Gla18$Gene %in% RPsGla$V1, ]
    GlaT24.RPs <- Gla24[Gla24$Gene %in% RPsGla$V1, ]
    #only common genes
    commonGla.RPs <- Reduce(x=list(GlaT12.RPs$Gene, GlaT18.RPs$Gene, GlaT24.RPs$Gene), f="intersect")
    GlaT12.RPs <- GlaT12.RPs[GlaT12.RPs$Gene %in% commonGla.RPs, ]
    GlaT18.RPs <- GlaT18.RPs[GlaT18.RPs$Gene %in% commonGla.RPs, ]
    GlaT24.RPs <- GlaT24.RPs[GlaT24.RPs$Gene %in% commonGla.RPs, ]
    
    
    ggplot() + xlim(c("Ecy12", "Ecy18", "Ecy24", "Eve12", "Eve18", "Eve24", "Gla12", "Gla18", "Gla24")) + 
      geom_boxplot(data=EcyT12.RPs, aes(x=1, y=log2FoldChange)) + 
      geom_boxplot(data=EcyT18.RPs, aes(x=2, y=log2FoldChange)) + 
      geom_boxplot(data=EcyT24.RPs, aes(x=3, y=log2FoldChange)) + 
      geom_boxplot(data=EveT12.RPs, aes(x=4, y=log2FoldChange)) + 
      geom_boxplot(data=EveT18.RPs, aes(x=5, y=log2FoldChange)) + 
      geom_boxplot(data=EveT24.RPs, aes(x=6, y=log2FoldChange)) + 
      geom_boxplot(data=GlaT12.RPs, aes(x=7, y=log2FoldChange)) + 
      geom_boxplot(data=GlaT18.RPs, aes(x=8, y=log2FoldChange)) + 
      geom_boxplot(data=GlaT24.RPs, aes(x=9, y=log2FoldChange)) + 
      theme_bw() + ggtitle("Ribosomal protein genes")
    
    ggsave("RPs.png")  
  
Ecy.RPs <- rbind(data.frame(logFC = EcyT12.RPs$log2FoldChange, temp = "T12"),
                 data.frame(logFC = EcyT18.RPs$log2FoldChange, temp = "T18"),
                 data.frame(logFC = EcyT24.RPs$log2FoldChange, temp = "T24"))
library(coin)
pairwise.wilcox.test(Ecy.RPs$logFC, Ecy.RPs$temp)


Eve.RPs <- rbind(data.frame(logFC = EveT12.RPs$log2FoldChange, temp = "T12"),
                 data.frame(logFC = EveT18.RPs$log2FoldChange, temp = "T18"),
                 data.frame(logFC = EveT24.RPs$log2FoldChange, temp = "T24"))
pairwise.wilcox.test(Eve.RPs$logFC, Eve.RPs$temp)

Gla.RPs <- rbind(data.frame(logFC = GlaT12.RPs$log2FoldChange, temp = "T12"),
                 data.frame(logFC = GlaT18.RPs$log2FoldChange, temp = "T18"),
                 data.frame(logFC = GlaT24.RPs$log2FoldChange, temp = "T24"))
pairwise.wilcox.test(Gla.RPs$logFC, Gla.RPs$temp)


## by DE genes
selectDEd <- function(tbl, threshold_p = 0.05, threshold_abslogFC = 1) {
  tbl[(tbl$log2FoldChange > threshold_abslogFC | tbl$log2FoldChange < -1 * threshold_abslogFC) & tbl$padj < threshold_p, ]
}

  EcyT12.DE <- selectDEd(Ecy12)
    EcyT12.DE[EcyT12.DE$Gene %in% intersect(EcyT12.DE$Gene, EcyT12.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EcyT12.DE$Gene, EcyT12.RPs$Gene))
  EcyT18.DE <- selectDEd(Ecy18); EcyT18.DE[EcyT18.DE$Gene %in% intersect(EcyT18.DE$Gene, EcyT18.RPs$Gene), ]
    EcyT18.DE[EcyT18.DE$Gene %in% intersect(EcyT18.DE$Gene, EcyT18.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EcyT18.DE$Gene, EcyT18.RPs$Gene))
  EcyT24.DE <- selectDEd(Ecy24); EcyT24.DE[EcyT24.DE$Gene %in% intersect(EcyT24.DE$Gene, EcyT24.RPs$Gene), ]
    EcyT24.DE[EcyT24.DE$Gene %in% intersect(EcyT24.DE$Gene, EcyT24.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EcyT24.DE$Gene, EcyT24.RPs$Gene))
  
    EveT12.DE <- selectDEd(Eve12)
    EveT12.DE[EveT12.DE$Gene %in% intersect(EveT12.DE$Gene, EveT12.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EveT12.DE$Gene, EveT12.RPs$Gene))
    EveT18.DE <- selectDEd(Eve18); EveT18.DE[EveT18.DE$Gene %in% intersect(EveT18.DE$Gene, EveT18.RPs$Gene), ]
    EveT18.DE[EveT18.DE$Gene %in% intersect(EveT18.DE$Gene, EveT18.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EveT18.DE$Gene, EveT18.RPs$Gene))
    EveT24.DE <- selectDEd(Eve24); EveT24.DE[EveT24.DE$Gene %in% intersect(EveT24.DE$Gene, EveT24.RPs$Gene), ]
    EveT24.DE[EveT24.DE$Gene %in% intersect(EveT24.DE$Gene, EveT24.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(EveT24.DE$Gene, EcyT24.RPs$Gene))
    
    
    
    GlaT12.DE <- selectDEd(Gla12)
    GlaT12.DE[GlaT12.DE$Gene %in% intersect(GlaT12.DE$Gene, GlaT12.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(GlaT12.DE$Gene, GlaT12.RPs$Gene))
    GlaT18.DE <- selectDEd(Gla18); GlaT18.DE[GlaT18.DE$Gene %in% intersect(GlaT18.DE$Gene, GlaT18.RPs$Gene), ]
    GlaT18.DE[GlaT18.DE$Gene %in% intersect(GlaT18.DE$Gene, GlaT18.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(GlaT18.DE$Gene, GlaT18.RPs$Gene))
    GlaT24.DE <- selectDEd(Gla24); GlaT24.DE[GlaT24.DE$Gene %in% intersect(GlaT24.DE$Gene, GlaT24.RPs$Gene), ]
    GlaT24.DE[GlaT24.DE$Gene %in% intersect(GlaT24.DE$Gene, GlaT24.RPs$Gene), c("Gene", "log2FoldChange", "padj")]
    length(intersect(GlaT24.DE$Gene, GlaT24.RPs$Gene))
    
    

# library(coin)
# wilcox.test(EcyT12.RPs$log2FoldChange, EcyT18.RPs$log2FoldChange)
# 
# wilcox.test(EveT12.RPs$log2FoldChange, EveT18.RPs$log2FoldChange)
# wilcox.test(EveT12.RPs$log2FoldChange, EveT24.RPs$log2FoldChange)
# 
# 
#     
#     HSPsEcy <- diamondEcy[grep("heat shock", diamondEcy$V14), ]
#     HSPsEve <- diamondEve[grep("heat shock", diamondEve$V14), ]
#     
#     EcyT12.hsps <- Ecy12[Ecy12$Gene %in% HSPsEcy$V1, ]
#     EcyT18.hsps <- Ecy18[Ecy18$Gene %in% HSPsEcy$V1, ]
#     EcyT24.hsps <- Ecy24[Ecy24$Gene %in% HSPsEcy$V1, ]
#     
#     EveT12.hsps <- Eve12[Eve12$Gene %in% HSPsEve$V1, ]
#     EveT18.hsps <- Eve18[Eve18$Gene %in% HSPsEve$V1, ]
#     EveT24.hsps <- Eve24[Eve24$Gene %in% HSPsEve$V1, ]
#     
#     ggplot() + xlim(c("Ecy12", "Ecy18", "Ecy24", "Eve12", "Eve18", "Eve24", "Gla12", "Gla18", "Gla24")) + 
#       geom_boxplot(data=EcyT12.hsps, aes(x=1, y=log2FoldChange)) + 
#       geom_boxplot(data=EcyT18.hsps, aes(x=2, y=log2FoldChange)) + 
#       geom_boxplot(data=EcyT24.hsps, aes(x=3, y=log2FoldChange)) +
#       geom_boxplot(data=EveT12.hsps, aes(x=4, y=log2FoldChange)) + 
#       geom_boxplot(data=EveT18.hsps, aes(x=5, y=log2FoldChange)) + 
#       geom_boxplot(data=EveT24.hsps, aes(x=6, y=log2FoldChange))
#     #   
#     
#     
#      Eve24$best.match.diamond <- sapply(as.character(Eve24$Gene), 
#                                function(x) diamondEve[diamondEve$V1 %in% x, "V14"])
#      write.csv(Eve24, "EveT24")
     
#     result$best.match.to.nr <- sapply(as.character(result$Gene), 
#                               function(x) fa[fa$Sequence.Name %in% x, "best.hit.to.nr"])
#     result$GO.Biological.Process <- sapply(as.character(result$Gene), 
#                       function(x) fa[fa$Sequence.Name %in% x, "GO.Biological.Process"])
#     result$GO.Molecular.Function <- sapply(as.character(result$Gene), 
#                      function(x) fa[fa$Sequence.Name %in% x, "GO.Molecular.Function"])
#     result$GO.Cellular.Component <- sapply(as.character(result$Gene), 
#                              function(x) fa[fa$Sequence.Name %in% x, "GO.Cellular.Component"])
#     
#     write.csv(result, paste0(comparison, ".", species, "common_DE.csv"))
#     return(result)
# 
# 
# 
# # #install.packages("openxlsx")
# # library(openxlsx)
# # write.xlsx(x = list("Ecy_ace_3h" = ecy_a03, "Ecy_ace_24h" = ecy_a24, "Ecy_phe_3h" = ecy_p03, "Ecy_phe_24h"= ecy_p24, 
# #                     "Ecy_ace+phe_3h" = ecy_ap03, "Ecy_ace+phe_24h" = ecy_ap24, 
# #                     "Eve_ace_24h" = eve_a24, "Eve_phe_3h" = eve_p03, "Eve_phe_24h"= eve_p24, 
# #                     "Eve_ace+phe_3h" = eve_ap03, "Eve_ace+phe_24h" = eve_ap24, 
# #                     "Gla_ace_3h" = gla_a03, "Gla_ace_24h" = gla_p24, "Gla_phe_3h" = gla_p03, "Gla_phe_24h"= gla_p24, 
# #                     "Gla_ace+phe_24h" = gla_ap24),
# #            file = "DE.xlsx")
# 
# 
# ## figure 1 (numer of DE genes)
# ## eg 
# ## table(ecy_ap24$salmon.edgeR.FC > 0)
# 
# 
# numberDE <- read.csv("./ace_phe_DE.csv")
# numberDE$samplegroup <- paste(numberDE$Species, numberDE$Time)
# 
# 
# ace <- numberDE[numberDE$Treatment=="ace", ]
# phe <- numberDE[numberDE$Treatment=="phe", ]
# ace_phe <- numberDE[numberDE$Treatment=="ace+phe", ]
# 
# 
# ggplot(ace, aes(x = X, fill=Species)) + 
#   geom_hline(yintercept = c(-80, -60, -40, -20, 20, 40, 60, 80), col = "grey") +
#   geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
#   geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
#   theme_classic(base_size = 16) + 
#   scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
#   ylab("Number of DE transcripts") + 
#   geom_hline(yintercept = 0, col = "grey") + 
#   ylim(-100, 100) 
# 
# ggsave("ace.svg",width=10, heigh=6, units = "cm")  
# 
# ggplot(phe, aes(x = X, fill=Species)) + 
#   geom_hline(yintercept = c(-140, -100, -60, -20, 20, 60, 100, 140, 180), col = "grey") +
#   geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
#   geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
#   theme_classic(base_size = 16) + 
#   scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
#  ylab("Number of DE transcripts") + 
#   geom_hline(yintercept = 0, col = "grey") + 
#   ylim(-150, 200) 
# 
# ggsave("phe.svg",width=10, heigh=6, units = "cm")  
# 
# 
# ggplot(ace_phe, aes(x = X, fill=Species)) + 
#   geom_hline(yintercept = c(-60, -40, -20, 20, 40, 60), col = "grey") +
#   geom_bar(aes(y=upregs), stat="identity", position="dodge") + 
#   geom_bar(aes(y=downregs), stat="identity", position="dodge") + 
#   theme_classic(base_size = 16) + 
#   scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
#   ylab("Number of DE transcripts") + 
#   geom_hline(yintercept = 0, col = "grey") + 
#   ylim(-75, 75) 
# 
# ggsave("ace_phe.svg",width=10, heigh=6, units = "cm")  