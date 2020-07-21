options(stringsAsFactors = F)
#library(UpSetR)
library(ggplot2)
library(reshape2)

## read files (made with salmon / DESeq2)
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

## read the annotation for each species
diamondEcy <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Ecy", "BCdTP1_ani.diamond.tsv"), head = F)
diamondEve <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Eve", "BCdTP1_ani.diamond.tsv"), head = F)
diamondGla <- read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Gla", "BCdTP1_ani.diamond.tsv"), head = F)

## find transcripts annotated as RPs
RPsEcy <- diamondEcy[grep("ribosomal pr", diamondEcy$V14), ]
RPsEve <- diamondEve[grep("ribosomal pr", diamondEve$V14), ]
RPsGla <- diamondGla[grep("ribosomal pr", diamondGla$V14), ]

## filter transcripts

EcyT12.RPs <- Ecy12[Ecy12$Gene %in% RPsEcy$V1, ]
EcyT18.RPs <- Ecy18[Ecy18$Gene %in% RPsEcy$V1, ]
EcyT24.RPs <- Ecy24[Ecy24$Gene %in% RPsEcy$V1, ]
#only common genes for the 3 lists
commonEcy.RPs <- Reduce(x=list(EcyT12.RPs$Gene, EcyT18.RPs$Gene, EcyT24.RPs$Gene), f="intersect")
EcyT12.RPs <- EcyT12.RPs[EcyT12.RPs$Gene %in% commonEcy.RPs, ]
EcyT18.RPs <- EcyT18.RPs[EcyT18.RPs$Gene %in% commonEcy.RPs, ]
EcyT24.RPs <- EcyT24.RPs[EcyT24.RPs$Gene %in% commonEcy.RPs, ]

EveT12.RPs <- Eve12[Eve12$Gene %in% RPsEve$V1, ]
EveT18.RPs <- Eve18[Eve18$Gene %in% RPsEve$V1, ]
EveT24.RPs <- Eve24[Eve24$Gene %in% RPsEve$V1, ]
#only common genes for the 3 lists
commonEve.RPs <- Reduce(x=list(EveT12.RPs$Gene, EveT18.RPs$Gene, EveT24.RPs$Gene), f="intersect")
EveT12.RPs <- EveT12.RPs[EveT12.RPs$Gene %in% commonEve.RPs, ]
EveT18.RPs <- EveT18.RPs[EveT18.RPs$Gene %in% commonEve.RPs, ]
EveT24.RPs <- EveT24.RPs[EveT24.RPs$Gene %in% commonEve.RPs, ]

GlaT12.RPs <- Gla12[Gla12$Gene %in% RPsGla$V1, ]
GlaT18.RPs <- Gla18[Gla18$Gene %in% RPsGla$V1, ]
GlaT24.RPs <- Gla24[Gla24$Gene %in% RPsGla$V1, ]
#only common genes for the 3 lists
commonGla.RPs <- Reduce(x=list(GlaT12.RPs$Gene, GlaT18.RPs$Gene, GlaT24.RPs$Gene), f="intersect")
GlaT12.RPs <- GlaT12.RPs[GlaT12.RPs$Gene %in% commonGla.RPs, ]
GlaT18.RPs <- GlaT18.RPs[GlaT18.RPs$Gene %in% commonGla.RPs, ]
GlaT24.RPs <- GlaT24.RPs[GlaT24.RPs$Gene %in% commonGla.RPs, ]

## plot (not used in the final version)
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

## make tables for each species 
Ecy.RPs <- rbind(data.frame(Species = "E. cyaneus", Gene = EcyT12.RPs$Gene, logFC = EcyT12.RPs$log2FoldChange, temperature = 12.4),
                 data.frame(Species = "E. cyaneus", Gene = EcyT18.RPs$Gene, logFC = EcyT18.RPs$log2FoldChange, temperature = 18.0),
                 data.frame(Species = "E. cyaneus", Gene = EcyT24.RPs$Gene, logFC = EcyT24.RPs$log2FoldChange, temperature = 24.4))

Eve.RPs <- rbind(data.frame(Species = "E. verrucosus", Gene = EveT12.RPs$Gene, logFC = EveT12.RPs$log2FoldChange, temperature = 12.4),
                 data.frame(Species = "E. verrucosus", Gene = EveT18.RPs$Gene, logFC = EveT18.RPs$log2FoldChange, temperature = 18.0),
                 data.frame(Species = "E. verrucosus", Gene = EveT24.RPs$Gene, logFC = EveT24.RPs$log2FoldChange, temperature = 24.4))

Gla.RPs <- rbind(data.frame(Species = "G. lacustris", Gene = GlaT12.RPs$Gene, logFC = GlaT12.RPs$log2FoldChange, temperature = 12.4),
                 data.frame(Species = "G. lacustris", Gene = GlaT18.RPs$Gene, logFC = GlaT18.RPs$log2FoldChange, temperature = 18.0),
                 data.frame(Species = "G. lacustris", Gene = GlaT24.RPs$Gene, logFC = GlaT24.RPs$log2FoldChange, temperature = 24.4))

## tests (not used in the final version)
library(coin)
pairwise.wilcox.test(Gla.RPs$logFC, Gla.RPs$temp)
pairwise.wilcox.test(Eve.RPs$logFC, Eve.RPs$temp)
pairwise.wilcox.test(Ecy.RPs$logFC, Ecy.RPs$temp)

## write RP data as in one multi-sheet xlsx file
library(openxlsx)
write.xlsx(list(Ecy.RPs, Eve.RPs, Gla.RPs), file = "RP_genes.xlsx")
