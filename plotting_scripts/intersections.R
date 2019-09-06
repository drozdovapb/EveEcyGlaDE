###Cadmium


{
  samples <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/samples_acute_stresses.csv", stringsAsFactors = F)
  
  samples <- samples[substr(samples$sample, 4, 5) == "Cd", ]
  
  numberDE <-data.frame(condition = samples$condition, species = samples$species, stringsAsFactors = F) 
  #remove replicates
  numberDE <- numberDE[!duplicated(numberDE), ]
  numberDE$upAll <- NA
  numberDE$downAll <- NA
  numberDE$upDecont <- NA
  numberDE$downDecont <- NA
  numberDE$upMetazoa <- NA
  numberDE$downMetazoa <- NA
  numberDE$samplegroup <- paste0(numberDE$species, numberDE$condition)
  
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    species <- thisrow[2]
    if (species == "Gla") { wdir <- paste0("/run/media/polina/Elements/transcriptome/DE/4-salmon/gla2") 
    } else {(wdir <- paste0("/run/media/polina/Elements/transcriptome/DE/4-salmon/", tolower(species))) }
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
    #numberDE$All[i] <- nrow(temp)
    #numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
    #numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
    #Ciliata excluded
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_decont_de.csv"))
    numberDE$Decont[i] <- nrow(temp)
    numberDE$upDecont[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downDecont[i] <- -sum(temp$log2FoldChange < 0)
    #Only Metazoa
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_metazoa_de.csv"))
    numberDE$Metazoa[i] <- nrow(temp)
    numberDE$upMetazoa[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downMetazoa[i] <- -sum(temp$log2FoldChange < 0)
  }
  
  
  library(ggplot2)
  
  str(numberDE)
  numberDE$samplegroup <- factor(numberDE$samplegroup, levels = rev(numberDE$samplegroup)) #mainly for plotting
  
  ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.5)) + 
    #geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
    #geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
    theme_classic(base_size = 16) + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
    geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
    ylab("number of DE genes") + 
    geom_hline(yintercept = 0, col = "white")
  
  
  ggsave(filename = "~/numberDE_Cd.svg", width = 125, height = 127, units = "mm")  
  #  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.svg")  
  #  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.pdf")  
}

{
  samples <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/samples_acute_stresses.csv", stringsAsFactors = F)
  
  samples <- samples[substr(samples$sample, 4, 5) == "Cd", ]
  
  numberDE <-data.frame(condition = samples$condition, species = samples$species, stringsAsFactors = F) 
  #remove replicates
  numberDE <- numberDE[!duplicated(numberDE), ]
  numberDE$upAll <- NA
  numberDE$downAll <- NA
  numberDE$upDecont <- NA
  numberDE$downDecont <- NA
  numberDE$upMetazoa <- NA
  numberDE$downMetazoa <- NA
  numberDE$samplegroup <- paste0(numberDE$species, numberDE$condition)
  
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    species <- thisrow[2]
    (wdir <- paste0("/run/media/polina/Elements/transcriptome/DE/9-reduced/", tolower(species))) 
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
    #numberDE$All[i] <- nrow(temp)
    #numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
    #numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
    #Ciliata excluded
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_decont_de.csv"))
    numberDE$Decont[i] <- nrow(temp)
    numberDE$upDecont[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downDecont[i] <- -sum(temp$log2FoldChange < 0)
    #Only Metazoa
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_metazoa_de.csv"))
    numberDE$Metazoa[i] <- nrow(temp)
    numberDE$upMetazoa[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downMetazoa[i] <- -sum(temp$log2FoldChange < 0)
  }
  
  
  
  str(numberDE)
  numberDE$samplegroup <- factor(numberDE$samplegroup, levels = rev(numberDE$samplegroup)) #mainly for plotting
  
  ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.5)) + 
    #geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
    #geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
    theme_classic(base_size = 16) + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
    geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
    ylab("number of DE genes") + 
    geom_hline(yintercept = 0, col = "white")
  
  
  ggsave(filename = "~/numberDE_Cd_red.svg", width = 125, height = 127, units = "mm")  
  #  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.svg")  
  #  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.pdf")  
}



options(stringsAsFactors = F)

library(UpSetR)
#eve and ecy are what we're interested at. And temperature and cadmium. 4-way is more or less okay
setwd("~/Documents/Paper1_stresses/data_tables/")
EveT24 <- read.csv("./reduced_each_by_own/EveLT1024_annot.csv")
EveT24up <- EveT24[(EveT24$log2FoldChange > 3 & EveT24$padj < 0.001), "best.nr.hit.diamond"]
writeLines(EveT24up, "~/EveT24up.txt")

EveT24down <- EveT24[(EveT24$log2FoldChange < -3 & EveT24$padj < 0.001), "best.nr.hit.diamond"]
writeLines(EveT24down, "~/EveT24down.txt")

EcyT24 <- read.csv("./reduced_each_by_own/EcyLT1024_annot.csv")
EcyT24up <- EcyT24[(EcyT24$log2FoldChange > 3 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]
EcyT24down <- EcyT24[(EcyT24$log2FoldChange < -3 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]

EveCd24 <- read.csv("./reduced_each_by_own/EveCd24_annot.csv")
EveCd24up <- EveCd24[EveCd24$log2FoldChange > 3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EveCd24up, "~/EveCd24up.txt")

EveCd24down <- EveCd24[EveCd24$log2FoldChange < -3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EveCd24down, "~/EveCd24down.txt")

EcyCd24 <- read.csv("./reduced_each_by_own/EcyCd24_annot.csv")
EcyCd24up <- EcyCd24[EcyCd24$log2FoldChange > 3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EcyCd24up, "~/EcyCd24up.txt")
EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange < -3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EcyCd24down, "~/EcyCd24down.txt")


intersect(EveT24up, EveCd24up)
intersect(EcyCd24up, EveCd24up)

intersect(intersect(EveT24up, EveCd24up), intersect(EcyCd24up, EveCd24up))
intersect(intersect(EveT24up, EveCd24up), intersect(EcyCd24up, EveCd24up))

Reduce(intersect, list(EveT24up, EcyT24up, EveCd24up, EcyCd24up))
#> Reduce(intersect, list(EveT24up, EcyT24up, EveCd24up, EcyCd24up))
#[1] "CAQ60114.1 70kDa heat shock protein [Gammarus locusta]"

EcyCd24[EcyCd24$best.nr.hit.diamond == "CAQ60114.1 70kDa heat shock protein [Gammarus locusta]", "gene"]
EcyT24[EcyT24$best.nr.hit.diamond == "CAQ60114.1 70kDa heat shock protein [Gammarus locusta]", "gene"]

intersect(EveT24down, EveCd24down)
intersect(EcyCd24down, EveCd24down)
intersect(intersect(EveT24down, EveCd24down), intersect(EcyCd24down, EveCd24down))


EveT24[EveT24$best.nr.hit.diamond == "AFI60316.1 heat shock protein 70 [Eulimnogammarus verrucosus]", "gene"]


GlaT24 <- read.csv("./reduced_each_by_own/GlaLT1024_annot.csv")
GlaT24up <- GlaT24[GlaT24$log2FoldChange > 3 & GlaT24$padj < 0.001, "best.nr.hit.diamond"]
GlaT24down <- GlaT24[GlaT24$log2FoldChange < -3 & GlaT24$padj < 0.001, "best.nr.hit.diamond"]

GlaCd24 <- read.csv("./reduced_each_by_own/GlaCd24_annot.csv")
GlaCd24up <- GlaCd24[GlaCd24$log2FoldChange > 3 & GlaCd24$padj < 0.001, "best.nr.hit.diamond"]
GlaCd24down <- GlaCd24[GlaCd24$log2FoldChange < -3 & GlaCd24$padj < 0.001, "best.nr.hit.diamond"]


intersect(EveT24up, GlaT24up)

writeLines(EcyCd24up, "~/EcyCd24up.txt")
EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange < -3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EcyCd24down, "~/EcyCd24down.txt")


intersect(EcyCd24up, EveCd24up)

intersect(EcyT24up, EveT24up)

#install.packages("Vennerable")
#library(Vennerable)
library(ggplot2)
USRinput <- list(EveLT10 = EveT24up, EcyLT10 = EcyT24up, EveCd = EveCd24up, EcyCd = EcyCd24up)
#what if we just 
svg()
upset(fromList(USRinput), keep.order = T, text.scale = 2) #, main.bar.color = "#009E73")
dev.off()




USRinput <- list(EveLT10 = EveT24up, EcyLT10 = EcyT24up, GlaLT10 = GlaT24up,
                 EveCd = EveCd24up, EcyCd = EcyCd24up, GlaCd = GlaCd24up)
#what if we just 
svg(width = 8, height = 6)
upset(fromList(USRinput), keep.order = T, text.scale = 2, nsets = 6)#, empty.intersections = T) #, main.bar.color = "#009E73")
dev.off()



intersect(GlaT24up, intersect(EcyT24up, EveT24up))





USRinput <- list(EveLT10 = EveT24down, EcyLT10 = EcyT24down, GlaLT10 = GlaT24down,
                 EveCd = EveCd24down, EcyCd = EcyCd24down, GlaCd = GlaCd24down)
svg(width = 8, height = 6)
upset(fromList(USRinput), keep.order = T, text.scale = 2, nsets = 6)
dev.off()

Reduce(intersect, list(EveT24down, EcyT24down, EveCd24down, EcyCd24down))
intersect(EveT24down, EcyT24down)
intersect(EveT24down, EveCd24down)

#library(venneuler)
install.packages("eulerr")
library(eulerr)
ve <- euler(c("Eve_10LT" = 28, "Eve_Cd" = 13, "Ecy_Cd" = 5, "Eve_10LT&Eve_Cd" = 3, 
          "Eve_10LT&Ecy_Cd = 1", "Eve_Cd&Ecy_Cd" = 1))

plot(ve)

# EveT24 <- read.csv("./ecy_reduced/EveT24_annot.csv")
# EveT24up <- EveT24[EveT24$log2FoldChange > 3 & EveT24$padj < 0.001, "best.nr.hit.diamond"]
# EveT24down <- EveT24[EveT24$log2FoldChange > 3 & EveT24$padj < 0.001, "best.nr.hit.diamond"]
# 
# EveCd24 <- read.csv("./ecy_reduced/EveCd24_annot.csv")
# EveCd24up <- EveCd24[EveCd24$log2FoldChange > 3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
# EveCd24down <- EveCd24[EveCd24$log2FoldChange > 3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
# 
# EcyCd24 <- read.csv("./ecy_reduced/EcyCd24_annot.csv")
# EcyCd24up <- EcyCd24[EcyCd24$log2FoldChange > 3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
# EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange > 3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
# 
# 
# intersect(EveT24up, EveCd24up)
# intersect(EcyCd24up, EveCd24up)




























library(UpSetR)
#eve and ecy are what we're interested at. And temperature and cadmium. 4-way is more or less okay
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced//")


extract_topN <- function(filename, up = T, ngenes = 50) {
  tbl <- read.csv(filename, stringsAsFactors = F)
  tbl_signif <- tbl[tbl$padj < 0.001, ]
  tbl_sorted <- tbl_signif[order(tbl_signif$log2FoldChange, decreasing = up) ,]
  tbl_top <- tbl_sorted[1:ngenes, ]
  minFC <- tbl_top[ngenes, "log2FoldChange"]
  message(paste0(minFC, " is the lowest FC in this list"))
  if (abs(minFC) < 1) on.exit("less than n genes")
  #return(tbl_top[,"best.hit.to.nr"])
  return(tbl_top[,"best.nr.hit.diamond"])
}

Eve10LT24down <- extract_topN("./EveLT1024_annot.csv", up = F, ngenes = 50)
EveCd24down <- extract_topN("./EveCd24_annot.csv", up = F, ngenes = 50)

intersect(Eve10LT24down, EveCd24down)
upset(fromList(list(VT = Eve10LT24down, VC = EveCd24down)))

EvePB24down <- extract_topN("./EvePB24_annot.csv", up = F, ngenes = 25)
upset(fromList(list(VT = Eve10LT24down, VA = EvePB24down)))

Ecy10LT24down <- extract_topN("./EcyLT1024_annot.csv", up = F, ngenes = 50)
upset(fromList(list(VT = Eve10LT24down, CT = Ecy10LT24down)))


EveT24up <- EveT24[EveT24$padj < 0.001, ][1:50,]
writeLines(EveT24up, "~/EveT24up.txt")

EveT24down <- EveT24[(EveT24$log2FoldChange < -3 & EveT24$padj < 0.001), "best.nr.hit.diamond"]
writeLines(EveT24down, "~/EveT24down.txt")

EcyT24 <- read.csv("./reduced_each_by_own/EcyLT1024_annot.csv")
EcyT24up <- EcyT24[(EcyT24$log2FoldChange > 3 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]
EcyT24down <- EcyT24[(EcyT24$log2FoldChange < -3 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]

EveCd24 <- read.csv("./reduced_each_by_own/EveCd24_annot.csv")
EveCd24up <- EveCd24[EveCd24$log2FoldChange > 3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EveCd24up, "~/EveCd24up.txt")

EveCd24down <- EveCd24[EveCd24$log2FoldChange < -3 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EveCd24down, "~/EveCd24down.txt")

EcyCd24 <- read.csv("./reduced_each_by_own/EcyCd24_annot.csv")
EcyCd24up <- EcyCd24[EcyCd24$log2FoldChange > 3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EcyCd24up, "~/EcyCd24up.txt")
EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange < -3 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
writeLines(EcyCd24down, "~/EcyCd24down.txt")



















library(UpSetR)
#eve and ecy are what we're interested at. And temperature and cadmium. 4-way is more or less okay
setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#setwd("~/Documents/Paper1_stresses/data_tables/ecy_reduced//")

extract_topN <- function(filename, up = T, ngenes = 50) {
  tbl <- read.csv(filename, stringsAsFactors = F)
  tbl_signif <- tbl[tbl$padj < 0.001, ]
  tbl_sorted <- tbl_signif[order(tbl_signif$log2FoldChange, decreasing = up) ,]
  tbl_top <- tbl_sorted[1:ngenes, ]
  minFC <- tbl_top[ngenes, "log2FoldChange"]
  message(paste0(minFC, " is the lowest FC in this list"))
  if (abs(minFC) < 1) on.exit("less than n genes")
  #return(tbl_top[,"best.hit.to.nr"])
  return(tbl_top[,"best.nr.hit.diamond"])
}

EvePB24down <- extract_topN("./EvePB24_annot.csv", up = F, ngenes = 10)
#EvePB24up <- extract_topN("./EvePB24_annot.csv", up = T, ngenes = 25)
EvePh24up <- extract_topN("./EvePh24_annot.csv", up = T, ngenes = 10)

upset(fromList(list(PBdown = EvePB24down, Phup = EvePh24up)))


setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#correlation?
EveAc <- read.csv("./EvePB24_annot.csv")
EvePh <- read.csv("./EvePh24_annot.csv")

merged <- merge(EveAc, EvePh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = merged$log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water") + ylab("Phenanthrene vs acetone") + 
  theme_gray(base_size = 16)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")
summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))

###Does it work with a different assembly? 
setwd("~/Documents/Paper1_stresses/data_tables/eve_deep/")
EveAc <- read.csv("./EvePB24_metazoa.csv")
EvePh <- read.csv("./EvePh24_metazoa.csv")
merged <- merge(EveAc, EvePh, by = "gene")
plot(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))

EveAc <- read.csv("./EvePB24_annot.csv")
EveT <- read.csv("./EveLT1024_annot.csv")
merged <- merge(EveAc, EveT, by = "gene")
plot(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
summary(fit)
plot(fit$residuals)
shapiro.test(sample(fit$residuals, 4999))
ks.test(fit$residuals, "pnorm")
res1=residuals(fit,type="response")
plot(res1)
ks.test(res1, "pnorm")


EveCd <- read.csv("./EveCd24_annot.csv")
merged <- merge(EveCd, EveT, by = "gene")
plot(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
summary(fit)



merged <- merge(EveT, EvePB, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = merged$log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
#  xlab("Acetone vs water") + ylab("Phenanthrene vs acetone") + 
  theme_gray(base_size = 16)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)


#Phenanthrene changes almost nothing
merged <- merge(EveT, EvePh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = merged$log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  #  xlab("Acetone vs water") + ylab("Phenanthrene vs acetone") + 
  theme_gray(base_size = 16)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)


merged <- merge(EvePh, EveAc, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = merged$log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  #  xlab("Acetone vs water") + ylab("Phenanthrene vs acetone") + 
  theme_gray(base_size = 16)
cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)








###Acetone and phenanthrene


setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#correlation?
EveAc <- read.csv("./EvePB24_annot.csv")
EvePh <- read.csv("./EvePh24_annot.csv")

merged <- merge(EveAc, EvePh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = merged$log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16)


ded <- merged[(abs(merged$log2FoldChange.x) > 1  & merged$padj.x < 0.001) |
              (abs(merged$log2FoldChange.y) > 1 & merged$padj.y < 0.001), ]
ggplot(data = ded, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. verrucosus, 24h")


cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")
cor(ded$log2FoldChange.x, ded$log2FoldChange.y, method = "pearson")
summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))

# ###Does it work with a different assembly? 
# setwd("~/Documents/Paper1_stresses/data_tables/eve_deep/")
# EveAc <- read.csv("./EvePB24_metazoa.csv")
# EvePh <- read.csv("./EvePh24_metazoa.csv")
# merged <- merge(EveAc, EvePh, by = "gene")
# plot(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
# summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))
# 
# EveAc <- read.csv("./EvePB24_annot.csv")
# EveT <- read.csv("./EveLT1024_annot.csv")
# merged <- merge(EveAc, EveT, by = "gene")
# plot(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
# cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "spearman")
# fit <- lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y)
# summary(fit)
# plot(fit$residuals)
# shapiro.test(sample(fit$residuals, 4999))
# ks.test(fit$residuals, "pnorm")
# res1=residuals(fit,type="response")
# plot(res1)
# ks.test(res1, "pnorm")


setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")
#correlation?
EcyAc <- read.csv("./EcyPB24_annot.csv")
EcyPh <- read.csv("./EcyPh24_annot.csv")

merged <- merge(EcyAc, EcyPh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 24h, DE only")


ded <- merged[(abs(merged$log2FoldChange.x) > 1  & merged$padj.x < 0.001) |
                (abs(merged$log2FoldChange.y) > 1 & merged$padj.y < 0.001), ]
ggplot(data = ded, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 24h, DE only")


cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")
cor(ded$log2FoldChange.x, ded$log2FoldChange.y, method = "pearson")
summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))


EveT24 <- read.csv("./EveLT1024_annot.csv")
EveT24up <- EcyT24[(EveT24$log2FoldChange > 1 & EveT24$padj < 0.001), "best.nr.hit.diamond"]
EveT24down <- EveT24[(EveT24$log2FoldChange < -1 & EveT24$padj < 0.001), "best.nr.hit.diamond"]

EveCd24 <- read.csv("./EveCd24_annot.csv")
EveCd24up <- EveCd24[EveCd24$log2FoldChange > 1 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
EveCd24down <- EveCd24[EveCd24$log2FoldChange < -1 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]

EveAcup <- EveAc[EveAc$log2FoldChange > 1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]
EveAcdown <- EveAc[EveAc$log2FoldChange < -1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]

USRinput <- list(EveT24up = EveT24up, EveCd24up = EveCd24up, EveAc24up = EveAcup)
upset(fromList(USRinput))

USRinput <- list(EveT24down = EveT24down, EveCd24down = EveCd24down, EveAc24down = EveAcdown)
upset(fromList(USRinput))






















####Other?

EcyAc <- read.csv("./EcyPB03_annot.csv")
EcyPh <- read.csv("./EcyPh03_annot.csv")

merged <- merge(EcyAc, EcyPh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 3h, DE only")

ded <- merged[(abs(merged$log2FoldChange.x) > 1  & merged$padj.x < 0.001) |
                (abs(merged$log2FoldChange.y) > 1 & merged$padj.y < 0.001), ]
ggplot(data = ded, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 3h, DE only")


cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson")
cor(ded$log2FoldChange.x, ded$log2FoldChange.y, method = "pearson")



#Correlation between T and Ac? 
EcyAc <- read.csv("./EcyPB24_annot.csv")
EcyT <- read.csv("./EcyLT1024_annot.csv")

merged <- merge(EcyAc, EcyT, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("Acetone vs water, log2FC") + ylab("LT10, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 3h, DE only")


#Correlation between T and Ph? 
EcyT <- read.csv("./EcyLT1024_annot.csv")
EcyPh <- read.csv("./EcyPh24_annot.csv")

merged <- merge(EcyT, EcyPh, by = "gene")
ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
  xlab("LT10, log2FC") + ylab("Phenanthrene vs acetone, log2FC") + 
  theme_gray(base_size = 16) + ggtitle("E. cyaneus, 3h, DE only")

####It's time to write up a function!

setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/")


getmecor <- function(xname, yname, x.lab, y.lab, plot.title) {
  xtable <- read.csv(paste0(xname, ".csv"))
  ytable <- read.csv(paste0(yname, ".csv"))
  merged <- merge(xtable, ytable, by = "gene")
  pm <- ggplot(data = merged, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
    xlab(x.lab) + ylab(y.lab) + 
    theme_gray(base_size = 16) + ggtitle(plot.title)
  print(pm)
  ggsave(paste0("../../multipanel_figures/acetone_phenanthrene/cor/", xname, "_vs_", yname, ".png"),
         width=5.77)
  
  ded <- merged[(abs(merged$log2FoldChange.x) > 1  & merged$padj.x < 0.001) |
                  (abs(merged$log2FoldChange.y) > 1 & merged$padj.y < 0.001), ]
  pd <- ggplot(data = ded, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
    geom_point(alpha = .3) + geom_smooth(method = "lm", col = "darkgray") + 
    xlab(x.lab) + ylab(y.lab) + 
    theme_gray(base_size = 16) + ggtitle(paste0(plot.title, ", DE only"))
  print(pd)
  
  print(cor(merged$log2FoldChange.x, merged$log2FoldChange.y, method = "pearson"))
  print(cor(ded$log2FoldChange.x, ded$log2FoldChange.y, method = "pearson"))
  #summary(lm(merged$log2FoldChange.x ~ merged$log2FoldChange.y))
}

getmecor("EcyPB24_annot", "EcyPh24_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "E. cyaneus, 24h")
getmecor("EcyPB03_annot", "EcyPh03_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "E. cyaneus, 3h")

getmecor("GlaPB24_annot", "GlaPh24_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "G. lacustris, 24h")
getmecor("GlaPB03_annot", "GlaPh03_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "G. lacustris, 3h")

getmecor("EvePB24_annot", "EvePh24_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "E. verrucosus, 24h")
getmecor("EvePB03_annot", "EvePh03_annot", 
         "Acetone vs water, log2FC", "Phenanthrene vs acetone, log2FC",
         "E. verrucosus, 3h")




getmecor("EvePB24_annot.csv", "EveLT1024_annot.csv", 
         "Acetone vs water, log2FC", "LT10, log2FC",
         "E. verrucosus, 24h")
getmecor("EvePB03_annot.csv", "EveLT1024_annot.csv", 
         "Acetone vs water, log2FC", "LT10, log2FC",
         "E. verrucosus, 24h")



getmecor("EveLT1024_annot.csv", "EvePh24_annot.csv", 
         "LT10, log2FC", "Phenanthrene vs acetone, log2FC",
         "E. verrucosus, 24h")



getmecor("GlaPB24_annot.csv", "GlaLT1024_annot.csv", 
         "Acetone vs water, log2FC", "LT10, log2FC",
         "G. lacustris, 24h")
getmecor("GlaLT1024_annot.csv", "GlaPh24_annot.csv", 
         "LT10, log2FC", "Phenanthrene vs acetone, log2FC",
         "G. lacustris, 24h")























EveT24 <- read.csv("./EveLT1024_annot.csv")
EveT24up <- EveT24[(EveT24$log2FoldChange > 1 & EveT24$padj < 0.001), "best.nr.hit.diamond"]
EveT24down <- EveT24[(EveT24$log2FoldChange < -1 & EveT24$padj < 0.001), "best.nr.hit.diamond"]

EveCd24 <- read.csv("./EveCd24_annot.csv")
EveCd24up <- EveCd24[EveCd24$log2FoldChange > 1 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]
EveCd24down <- EveCd24[EveCd24$log2FoldChange < -1 & EveCd24$padj < 0.001, "best.nr.hit.diamond"]

EveAc <- read.csv("./EvePB24_annot.csv")
EveAcup <- EveAc[EveAc$log2FoldChange > 1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]
EveAcdown <- EveAc[EveAc$log2FoldChange < -1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]

intersect(intersect(EveT24down, EveCd24down), EveAcdown)
intersect(intersect(EveT24up, EveCd24up), EveAcup)


EcyT24 <- read.csv("./EcyLT1024_annot.csv")
EcyT24up <- EcyT24[(EcyT24$log2FoldChange > 1 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]
EcyT24down <- EcyT24[(EcyT24$log2FoldChange < -1 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]

EcyCd24 <- read.csv("./EcyCd24_annot.csv")
EcyCd24up <- EcyCd24[EcyCd24$log2FoldChange > 1 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange < -1 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]

EcyAc <- read.csv("./EcyPB24_annot.csv")
EcyAcup <- EcyAc[EcyAc$log2FoldChange > 1 & EcyAc$padj < 0.001, "best.nr.hit.diamond"]
EcyAcdown <- EcyAc[EcyAc$log2FoldChange < -1 & EcyAc$padj < 0.001, "best.nr.hit.diamond"]

intersect(EcyT24down, EcyCd24down)

intersect(intersect(EcyT24down, EcyCd24down), EcyAcdown)
intersect(intersect(EcyT24up, EcyCd24up), EcyAcup)


















EveT24 <- read.csv("./EveLT1024_annot.csv")
EveT24up <- EveT24[(EveT24$log2FoldChange > 1 & EveT24$padj < 0.001), "gene"]
EveT24down <- EveT24[(EveT24$log2FoldChange < -1 & EveT24$padj < 0.001), "gene"]

EveCd24 <- read.csv("./EveCd24_annot.csv")
EveCd24up <- EveCd24[EveCd24$log2FoldChange > 1 & EveCd24$padj < 0.001, "gene"]
EveCd24down <- EveCd24[EveCd24$log2FoldChange < -1 & EveCd24$padj < 0.001, "gene"]

upset(data = fromList(list(tu = EveT24up, td = EveT24down, cu = EveCd24up, cd = EveCd24down)))


EveAc <- read.csv("./EvePB24_annot.csv")
EveAcup <- EveAc[EveAc$log2FoldChange > 1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]
EveAcdown <- EveAc[EveAc$log2FoldChange < -1 & EveAc$padj < 0.001, "best.nr.hit.diamond"]

intersect(intersect(EveT24down, EveCd24down), EveAcdown)
intersect(intersect(EveT24up, EveCd24up), EveAcup)


EcyT24 <- read.csv("./EcyLT1024_annot.csv")
EcyT24up <- EcyT24[(EcyT24$log2FoldChange > 1 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]
EcyT24down <- EcyT24[(EcyT24$log2FoldChange < -1 & EcyT24$padj < 0.001), "best.nr.hit.diamond"]

EcyCd24 <- read.csv("./EcyCd24_annot.csv")
EcyCd24up <- EcyCd24[EcyCd24$log2FoldChange > 1 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]
EcyCd24down <- EcyCd24[EcyCd24$log2FoldChange < -1 & EcyCd24$padj < 0.001, "best.nr.hit.diamond"]

EcyAc <- read.csv("./EcyPB24_annot.csv")
EcyAcup <- EcyAc[EcyAc$log2FoldChange > 1 & EcyAc$padj < 0.001, "best.nr.hit.diamond"]
EcyAcdown <- EcyAc[EcyAc$log2FoldChange < -1 & EcyAc$padj < 0.001, "best.nr.hit.diamond"]