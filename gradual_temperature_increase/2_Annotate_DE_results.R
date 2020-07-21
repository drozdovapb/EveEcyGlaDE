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

## Ecy
Ecy12$best.match.diamond <- sapply(as.character(Ecy12$Gene), function(x) diamondEcy[diamondEcy$V1 %in% x, "V14"])
Ecy18$best.match.diamond <- sapply(as.character(Ecy18$Gene), function(x) diamondEcy[diamondEcy$V1 %in% x, "V14"])
Ecy24$best.match.diamond <- sapply(as.character(Ecy24$Gene), function(x) diamondEcy[diamondEcy$V1 %in% x, "V14"])
## Eve
Eve12$best.match.diamond <- sapply(as.character(Eve12$Gene), function(x) diamondEve[diamondEve$V1 %in% x, "V14"])
Eve18$best.match.diamond <- sapply(as.character(Eve18$Gene), function(x) diamondEve[diamondEve$V1 %in% x, "V14"])
Eve24$best.match.diamond <- sapply(as.character(Eve24$Gene), function(x) diamondEve[diamondEve$V1 %in% x, "V14"])
## Gla
Gla12$best.match.diamond <- sapply(as.character(Gla12$Gene), function(x) diamondGla[diamondGla$V1 %in% x, "V14"])
Gla18$best.match.diamond <- sapply(as.character(Gla18$Gene), function(x) diamondGla[diamondGla$V1 %in% x, "V14"])
Gla24$best.match.diamond <- sapply(as.character(Gla24$Gene), function(x) diamondGla[diamondGla$V1 %in% x, "V14"])


write.csv(Ecy12, "Ecy.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results.annot")
write.csv(Ecy18, "Ecy.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results.annot")
write.csv(Ecy24, "Ecy.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results.annot")

write.csv(Eve12, "Eve.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results.annot")
write.csv(Eve18, "Eve.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results.annot")
write.csv(Eve24, "Eve.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results.annot")

write.csv(Gla12, "Gla.isoform.counts.matrix.T12_vs_B6.DESeq2.DE_results.annot")
write.csv(Gla18, "Gla.isoform.counts.matrix.T18_vs_B6.DESeq2.DE_results.annot")
write.csv(Gla24, "Gla.isoform.counts.matrix.T24_vs_B6.DESeq2.DE_results.annot")
