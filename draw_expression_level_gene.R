options(stringsAsFactors = F)

#read data.
#for now verrucosus
dat <- read.delim("/homes/brauerei/polina/Documents/Paper2_time_series/eve_all.cnt", sep = " ")
dat <- dat[,grep("deep", names(dat), invert = T)]


#grep(dat, "")

library(ggplot2)
library(reshape2)

datmelt <- melt(dat, id.vars = "Name", 
                variable.name = "Sample", value.name = "TPM", factorsAsStrings = T)






diamondv <- read.delim("~/Documents/transcriptome_annotation/EveBCdTP1_ani.diamond.tsv", head = F,
                       stringsAsFactors = F)

thesecontigs <- diamondv[grep("heat shock 70", diamondv$V14), c("V1", "V14")]
hsp70 <- grep("heat shock", diamondv$V14) & grep("70", diamondv$V14)
thesecontigs <- diamondv[grep("heat shock 70", diamondv$V14), c("V1", "V14")] 





drawmeplease <- function(term) {
      
  thesecontigs <- diamondv[grep(term, diamondv$V14), c("V1", "V14")] 
  print(thesecontigs)
  
    for (i in 1:nrow(thesecontigs)) {
    thisdat <- dat[dat$Name %in% thesecontigs[i,1], ]
    thisdat$Description <- thesecontigs[i,2]
    datmelt <- melt(thisdat, id.vars = c("Name", "Description"), 
                    variable.name = "Sample", value.name = "TPM", factorsAsStrings = T)
    
    
    #intersect(dat$Name, diamondv$V1)
    
    datmelt$sample <- as.character(datmelt$Sample)
    datmelt$samplegroup <- sapply(datmelt$sample, function(x) strsplit(x, "_")[[1]][2])
    datmelt$condition <- substr(datmelt$samplegroup, 4, 20)
    datmelt$fcondition <- factor(datmelt$condition, levels = c("B6", "B12", "B18", "B24", "B3", 
                                                               "T12", "T18", "T24", "10LT3", "10LT24",
                                                               "Cd3", "Cd24", "PB3", "PB24", "P3", "P24"))
    
    p <- ggplot(data = datmelt, aes(x = fcondition, y = TPM)) + 
      geom_boxplot() + geom_point() + ggtitle(paste(datmelt$Name[1], "\n", datmelt$Description[1]))
    print(p)  
  }
}
  

drawmeplease("70kDa heat shock")
drawmeplease("catalase")

drawmeplease("peroxidase")
drawmeplease("thioredoxin")

#glutathione S-transferase; peroxidasin; selenoprotein W; thioredoxin; superoxide dismutase. 

drawmeplease("glutathione S-transferase")


drawmeplease("vitellogenin")
