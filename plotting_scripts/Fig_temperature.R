##Temperature

##LT10
##Basically, this just plots GraphPad values from Till
{library(ggplot2)
  library(reshape2)

  pred <- read.csv("../../temperature/predictions.csv")
  names(pred) <- c("Temperature", "Eve", "Ecy", "Gla")
  
  dat <- read.delim("../../temperature/temp-mort.csv")
  dat$Temperature <- as.numeric(dat$Temperature)
  
pLT10 <-  ggplot(dat, aes(x=Temperature, Mortality)) + 
    #geom_point(aes(col = Species)) +
    geom_count(aes(col = Species)) + 
    scale_radius(range=c(3,6), breaks=1:3, name = "Observations") + #scale_size_area(max_size = 6)
    scale_color_manual(values = c("#56B4E9", "#009E73", "#E69F00"),
                       labels = c("E. cyaneus", "E. verrucosus", "G. lacustris")) + #, 
    # labels = c(expression(italic("E. verrucosus")),
    #            expression(italic("E. cyaneus")), 
    #            expression(italic("G. lacustris")))) +
    theme_classic(base_size = 16)  + 
    geom_line(data = pred, aes(x=Temperature, y = Eve), col = "#009E73") +
    geom_line(data = pred, aes(x=Temperature, y = Ecy), col = "#56B4E9") +
    geom_line(data = pred, aes(x=Temperature, y = Gla), col = "#E69F00") + 
    geom_hline(aes(yintercept = 10, linetype = "LT10")) + 
    geom_hline(aes(yintercept = 50, linetype = "LT50")) +
    scale_linetype_manual(name = "", values = c("dotted", "dotdash"), 
                          guide = guide_legend(override.aes = list(lty = c("dotted", "dotdash")))) + 
    scale_x_continuous(breaks = seq(20, 29, 1)) + 
    scale_y_continuous(breaks = c(0, 10, 25, 50, 75, 100)) + 
    guides(color = guide_legend(order = -1, 
                                label.theme = element_text(face = "italic", size = 14),
                                title.theme = element_text(size = 16)),
           size = guide_legend(size = 1, title.theme = element_text(size = 16),
                               label.theme = element_text(size = 12))) + 
    xlab("Temperature, °C") + ylab("Mortality, %")
}

pLT10
#ggsave(filename = "./LT10_50.svg", width = 7)
ggsave(filename = "./LT10_50.svg", width = 8)














#Numbers of DE genes
{
  samples <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/samples_acute_stresses.csv", stringsAsFactors = F)
  
  samples <- samples[substr(samples$sample, 4, 7) == "10LT", ]
  
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
  
  
  ggsave(filename = "~/numberDE_temp.svg", width = 125, units = "mm")  
#  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.svg")  
#  ggsave(filename = "~/Documents/Paper1_stresses/temperature/numberDE.pdf")  
  }








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



#Funny... if all by Eve?
{
  for (i in 1:nrow(numberDE)) {
  thisrow <- numberDE[i, ]
  wdir <- "/run/media/polina/Elements/transcriptome/DE/4-salmon/eve"
  setwd(wdir)
  #All genes (no filtering)
  temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
  numberDE$All[i] <- nrow(temp)
  numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
  numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
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

ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.33)) + 
  geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
  geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
  theme_light() + 
  coord_flip() + 
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
  geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
  geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
  geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
  geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
  ylab("number of DE genes") + 
  geom_hline(yintercept = 0, col = "white") 
}
#Not the same but similar

#Eve_deep?
{
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    wdir <- "/run/media/polina/Elements/transcriptome/DE/4-salmon/eve_deep"
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
    numberDE$All[i] <- nrow(temp)
    numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
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
  
  ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.33)) + 
    geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
    geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
    theme_light() + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
    geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
    ylab("number of DE genes") + 
    geom_hline(yintercept = 0, col = "white") 
}
#Something funny happens... anyway, the same

#Ecy?
{
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    wdir <- "/run/media/polina/Elements/transcriptome/DE/4-salmon/ecy"
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
    numberDE$All[i] <- nrow(temp)
    numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
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
  
  ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.33)) + 
    geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
    geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
    theme_light() + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
    geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
    ylab("number of DE genes") + 
    geom_hline(yintercept = 0, col = "white") 
}
#Obviously the same

#Gla?
{
  for (i in 1:nrow(numberDE)) {
    thisrow <- numberDE[i, ]
    wdir <- "/run/media/polina/Elements/transcriptome/DE/4-salmon/gla2"
    setwd(wdir)
    #All genes (no filtering)
    temp <- read.csv(paste0(thisrow[2], thisrow[1], "_annot_de.csv"))
    numberDE$All[i] <- nrow(temp)
    numberDE$upAll[i] <- sum(temp$log2FoldChange > 0)
    numberDE$downAll[i] <- -sum(temp$log2FoldChange < 0)
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
  
  ggplot(numberDE, aes(x = samplegroup, fill=species, alpha = 0.33)) + 
    geom_bar(aes(y=upAll), stat="identity", position="dodge", col = "black") + 
    geom_bar(aes(y=downAll), stat="identity", position="dodge", col = "black") + 
    theme_light() + 
    coord_flip() + 
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) + 
    geom_bar(aes(y=upDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=downDecont), stat="identity", position="dodge") + 
    geom_bar(aes(y=upMetazoa), stat="identity", position="dodge") + 
    geom_bar(aes(y=downMetazoa), stat="identity", position="dodge") + 
    ylab("number of DE genes") + 
    geom_hline(yintercept = 0, col = "white") 
}
#The same, even though the scale is really different and response of Ecy and Gla is the same!!!!