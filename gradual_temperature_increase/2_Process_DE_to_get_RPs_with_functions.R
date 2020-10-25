options(stringsAsFactors = F)
#library(UpSetR)
library(ggplot2)
library(reshape2)


## read the annotation for each species (for the sake of time)
diamond <- list (Ecy = read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Ecy", "BCdTP1_ani.diamond.tsv"), head = F),
                 Eve = read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Eve", "BCdTP1_ani.diamond.tsv"), head = F),
                 Gla = read.delim(paste0("~/Research/Projects/DE/2-Assemblies_and_annotation/annotation/",
                                "Gla", "BCdTP1_ani.diamond.tsv"), head = F))

## First, filter the RP-encoding transcripts from the whole DE table
get_RPs <- function(Species, Condition){
  ## read files (made with salmon / DESeq2)
  de.result <- read.delim(paste0(Species, ".isoform.counts.matrix.", Condition, "_vs_B6.DESeq2.DE_results"))
  de.result$Gene <- row.names(de.result)
  ## find transcripts annotated as RPs
  RPs <- diamond[[Species]][grep("ribosomal pr", diamond[[Species]]$V14), ]
  ## filter transcripts
  result.RPs <- de.result[de.result$Gene %in% RPs$V1, ]
  return(result.RPs)
}

## Get only common RPs for each species
get_common_RPs <- function(species) {
  T12.RPs <- get_RPs(species, "T12")
  T18.RPs <- get_RPs(species, "T18")
  T24.RPs <- get_RPs(species, "T24")
  #only common genes for the 3 lists
  common.RPs <- Reduce(x=list(T12.RPs$Gene, T18.RPs$Gene, T24.RPs$Gene), f="intersect")
  T12.RPs <- T12.RPs[T12.RPs$Gene %in% common.RPs, ]
  T18.RPs <- T18.RPs[T18.RPs$Gene %in% common.RPs, ]
  T24.RPs <- T24.RPs[T24.RPs$Gene %in% common.RPs, ]
 common.RPs <- rbind(data.frame(Species = species, Gene = T12.RPs$Gene, logFC = T12.RPs$log2FoldChange, temperature = 12.4),
                  data.frame(Species = species, Gene = T18.RPs$Gene, logFC = T18.RPs$log2FoldChange, temperature = 18.0),
                  data.frame(Species = species, Gene = T24.RPs$Gene, logFC = T24.RPs$log2FoldChange, temperature = 24.4))
 return(common.RPs)
}

Ecy.RPs <- get_common_RPs("Ecy")
Eve.RPs <- get_common_RPs("Eve")
Gla.RPs <- get_common_RPs("Gla")

## plot (not used in the final version)
ggplot() + xlim(c("Ecy12", "Ecy18", "Ecy24", "Eve12", "Eve18", "Eve24", "Gla12", "Gla18", "Gla24")) + 
  geom_boxplot(data=Ecy.RPs[Ecy.RPs$temperature==12.4,], aes(x=1, y=logFC)) + 
  geom_boxplot(data=Ecy.RPs[Ecy.RPs$temperature==18.0,], aes(x=2, y=logFC)) + 
  geom_boxplot(data=Ecy.RPs[Ecy.RPs$temperature==24.4,], aes(x=3, y=logFC)) + 
  geom_boxplot(data=Eve.RPs[Eve.RPs$temperature==12.4,], aes(x=4, y=logFC)) + 
  geom_boxplot(data=Eve.RPs[Eve.RPs$temperature==18.0,], aes(x=5, y=logFC)) + 
  geom_boxplot(data=Eve.RPs[Eve.RPs$temperature==24.4,], aes(x=6, y=logFC)) + 
  geom_boxplot(data=Gla.RPs[Gla.RPs$temperature==12.4,], aes(x=7, y=logFC)) + 
  geom_boxplot(data=Gla.RPs[Gla.RPs$temperature==18.0,], aes(x=8, y=logFC)) + 
  geom_boxplot(data=Gla.RPs[Gla.RPs$temperature==24.4,], aes(x=9, y=logFC)) + 
  theme_bw() + ggtitle("Ribosomal protein genes")
ggsave("RPs.png")  

## tests (not used in the final version)
library(coin)
pairwise.wilcox.test(Gla.RPs$logFC, Gla.RPs$temp)
pairwise.wilcox.test(Eve.RPs$logFC, Eve.RPs$temp)
pairwise.wilcox.test(Ecy.RPs$logFC, Ecy.RPs$temp)

## write RP data as in one multi-sheet xlsx file
library(openxlsx)
write.xlsx(list(Ecy.RPs, Eve.RPs, Gla.RPs), file = "RP_genes.xlsx")
