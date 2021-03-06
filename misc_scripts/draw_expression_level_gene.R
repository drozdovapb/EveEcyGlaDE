options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)

#####
## this part is to be run only once to filter the sets
setwd("/homes/brauerei/polina/Documents/Paper2_time_series/")
ts_samples <- c("B06", "B12", "B18", "B24", "T12", "T18", "T24")
samples <- read.csv("samples.csv")
#Eve
eve <- read.delim("./eve_all.cnt", sep = " ")
eve_ts_sp <- samples[samples$species == "Eve" & samples$condition %in% ts_samples, "sample"]
evets <- eve[, c(1, which(names(eve) %in% eve_ts_sp))]
write.csv(evets, "eve_T.cnt.csv", row.names = F)
#Ecy
ecy <- read.delim("./ecy_all.cnt", sep = " ")
ecy_ts_sp <- samples[samples$species == "Ecy" & samples$condition %in% ts_samples, "sample"]
ecyts <- ecy[, c(1, which(names(ecy) %in% ecy_ts_sp))]
write.csv(ecyts, "ecy_T.cnt.csv", row.names = F)
#Gla
gla <- read.delim("./gla_all.cnt", sep = " ")
gla_ts_sp <- samples[samples$species == "Gla" & samples$condition %in% ts_samples, "sample"]
glats <- gla[, c(1, which(names(gla) %in% gla_ts_sp))]
write.csv(glats, "gla_T.cnt.csv", row.names = F)


## the same about this part, for acute stressors
setwd("/homes/brauerei/polina/Documents/Paper1_stresses/")
paper1_samples <- c("B03", "B24h", "LT1003", "LT1024", "Cd03", "Cd24")
samples <- read.csv("samples.csv")
#Eve
eve <- read.delim("../Paper2_time_series/eve_all.cnt", sep = " ")
eve_ac_sp <- samples[samples$species == "Eve" & samples$condition %in% paper1_samples, "sample"]
evea <- eve[, c(1, which(names(eve) %in% eve_ac_sp))]
write.csv(evea, "eve_paper1.cnt.csv", row.names = F)
#Ecy
ecy <- read.delim("../Paper2_time_series/ecy_all.cnt", sep = " ")
ecy_ac_sp <- samples[samples$species == "Ecy" & samples$condition %in% paper1_samples, "sample"]
ecya <- ecy[, c(1, which(names(ecy) %in% ecy_ac_sp))]
write.csv(ecya, "ecy_paper1.cnt.csv", row.names = F)
#Gla
gla <- read.delim("../Paper2_time_series/gla_all.cnt", sep = " ")
gla_ac_sp <- samples[samples$species == "Gla" & samples$condition %in% paper1_samples, "sample"]
glaa <- gla[, c(1, which(names(gla) %in% gla_ac_sp))]
write.csv(glaa, "gla_paper1.cnt.csv", row.names = F)



#####

#read data.
#for now verrucosus


#if I use the verrucosus data and want to get rid of the deep sampels...
#dat <- dat[,grep("deep", names(dat), invert = T)]


## a function to draw by term (if you don't know the particular contig)
draw_by_term <- function(term, dat = datv, annotation = diamondv, species = "Eve") {
  ## get part of the data for the particular term
  thesecontigs <- annotation[grep(term, annotation$V14), c("V1", "V14")] 
  print(thesecontigs)  # to get an idea of what we have
  #intersect(dat$Name, diamondv$V1) #sanity check 
  ## and now 
  for (i in 1:nrow(thesecontigs)) {
    thisdat <- dat[dat$Name %in% thesecontigs[i,1], ]
    thisdat$Description <- thesecontigs[i,2]
    datmelt <- melt(thisdat, id.vars = c("Name", "Description"), 
                    variable.name = "Sample", value.name = "TPM", factorsAsStrings = T)
    datmelt$sample <- as.character(datmelt$Sample)
    datmelt$samplegroup <- sapply(datmelt$sample, function(x) strsplit(x, "_")[[1]][1])
    datmelt$condition <- substr(datmelt$samplegroup, 4, 20)
    datmelt$fcondition <- factor(datmelt$condition, levels = c("B3", "B6", "B12", "B18", "B24", 
                                                               "T12", "T18", "T24", "10LT3", "10LT24",
                                                               "Cd3", "Cd24", "PB3", "PB24", "P3", "P24"))
    #datmelt$fcondition <- factor(datmelt$condition, levels = c("B6", "B12", "B18", "B24",  
    #                                                           "T12", "T18", "T24"))
    title <- paste(species, datmelt$Name[1], "\n", datmelt$Description[1])
    if (grepl("/", datmelt$Description[1])) {
      safe_description <- gsub("/", "-", datmelt$Description[1])
    } else safe_description <- datmelt$Description[1]
    fname <- paste0(species, "_", datmelt$Name[1], "_", safe_description, ".png")
    p <- ggplot(data = datmelt, aes(x = fcondition, y = TPM)) + 
      geom_boxplot() + geom_point() + ggtitle(title)
    print(p)
    ggsave(paste0(fname), width = 20, height = 15, units = "cm")
  }
}

## and a function to draw by id if you know it already
draw_by_id <- function(id, dat = datv, annotation = diamondv, species = "Eve") {
  ## get part of the data for the particular term
  thesecontigs <- annotation[grep(id, annotation$V1), c("V1", "V14")] 
  print(thesecontigs)  # to get an idea of what we have
  #intersect(dat$Name, diamondv$V1) #sanity check 
  ## and now 
  for (i in 1:nrow(thesecontigs)) {
    thisdat <- dat[dat$Name %in% thesecontigs[i,1], ]
    thisdat$Description <- thesecontigs[i,2]
    datmelt <- melt(thisdat, id.vars = c("Name", "Description"), 
                    variable.name = "Sample", value.name = "TPM", factorsAsStrings = T)
    datmelt$sample <- as.character(datmelt$Sample)
    datmelt$samplegroup <- sapply(datmelt$sample, function(x) strsplit(x, "_")[[1]][1])
    datmelt$condition <- substr(datmelt$samplegroup, 4, 20)
    datmelt$fcondition <- factor(datmelt$condition, levels = c("B3", "B6", "B12", "B18", "B24", 
                                                               "T12", "T18", "T24", "10LT3", "10LT24",
                                                               "Cd3", "Cd24", "PB3", "PB24", "P3", "P24"))
    #datmelt$fcondition <- factor(datmelt$condition, levels = c("B6", "B12", "B18", "B24", 
    #                                                           "T12", "T18", "T24"))
    title <- paste(species, datmelt$Name[1], "\n", datmelt$Description[1])
    ## prevent problems with filenames
    if (grepl("/", datmelt$Description[1])) {
      safe_description <- gsub("/", "-", datmelt$Description[1])
    } else safe_description <- datmelt$Description[1]
    fname <- paste0(species, "_", datmelt$Name[1], "_", safe_description, ".png")
    p <- ggplot(data = datmelt, aes(x = fcondition, y = TPM)) + 
      geom_boxplot() + geom_point() + ggtitle(title)
    print(p)
    ggsave(paste0(fname), width = 20, height = 15, units = "cm")
  }
}



## Some examples
drawmeplease("70kDa heat shock")
drawmeplease("catalase")
drawmeplease("peroxidase")
drawmeplease("thioredoxin")
#glutathione S-transferase; peroxidasin; selenoprotein W; thioredoxin; superoxide dismutase. 
drawmeplease("glutathione S-transferase")
drawmeplease("vitellogenin")
drawmeplease("cytochrome c oxidase subunit I")
drawmeplease("citrate")


setwd("/homes/brauerei/polina/Documents/Paper2_time_series/")

## But so far we're the most interested in gradual temperature increase
#evedat <- read.delim("./eve_all.cnt", sep = " ")
evedat <- read.csv("./eve_T.cnt.csv")
ecydat <- read.csv("./ecy_T.cnt.csv")
gladat <- read.csv("./gla_T.cnt.csv")




eveALLdat <- read.csv("./eve_all.cnt", sep = " ")
ecyALLdat <- read.csv("./ecy_all.cnt", sep = " ")
glaALLdat <- read.csv("./gla_all.cnt", sep = " ")

#paper 1 data
setwd("/homes/brauerei/polina/Documents/Paper1_stresses/counts/")
eveTCdat <- read.csv("./eve_paper1.cnt.csv")
ecyTCdat <- read.csv("./ecy_paper1.cnt.csv")
glaTCdat <- read.csv("./gla_paper1.cnt.csv")

diamondv <- read.delim("~/Documents/transcriptome_annotation/EveBCdTP1_ani.diamond.tsv", head = F,
                       stringsAsFactors = F)
diamondc <- read.delim("~/Documents/transcriptome_annotation/EcyBCdTP1_ani.diamond.tsv", head = F,
                       stringsAsFactors = F)
diamondl <- read.delim("~/Documents/transcriptome_annotation/GlaBCdTP1_ani.diamond.tsv", head = F,
                       stringsAsFactors = F)

drawmeplease("heat")


draw_by_id("TRINITY_DN371471", dat = evedat, annotation=diamondv, species = "Eve")
draw_by_id("TRINITY_DN510477", dat = ecydat, annotation=diamondc, species = "Ecy")
draw_by_id("TRINITY_DN365602", dat = gladat, annotation=diamondl, species = "Gla")


draw_by_term("cytochrome c oxidase", dat = evedat, annotation=diamondv, species = "Eve")
draw_by_term("cytochrome c oxidase", dat = ecydat, annotation=diamondc, species = "Ecy")
draw_by_term("cytochrome c oxidase", dat = gladat, annotation=diamondl, species = "Gla")

#citrate synthase
draw_by_term("citrate synthase", dat = evedat, annotation=diamondv, species = "Eve")
draw_by_term("citrate synthase", dat = ecydat, annotation=diamondc, species = "Ecy")
draw_by_term("citrate synthase", dat = gladat, annotation=diamondl, species = "Gla")


draw_by_term("pyruvate kinase", dat = evedat, annotation=diamondv, species = "Eve")



draw_by_term("vitellogenin", dat = evedat, annotation=diamondv, species = "Eve")
draw_by_term("vitellogenin", dat = ecydat, annotation=diamondc, species = "Ecy")
draw_by_term("vitellogenin", dat = gladat, annotation=diamondl, species = "Gla")


#COX3 in Eve
draw_by_id("TRINITY_DN280825_c0_g1_i1", dat = eveALLdat, annotation = diamondv, species = "Eve")


draw_by_id("TRINITY_DN280825_c0_g1_i1", dat = eveTCdat, annotation = diamondv, species = "Eve")

draw_by_id("TRINITY_DN280825", dat = eveTCdat, annotation = diamondv, species = "Eve")

draw_by_term("COX3", dat = eveTCdat, annotation = diamondv, species = "Eve")
draw_by_term("COX3", dat = ecyTCdat, annotation = diamondc, species = "Ecy")
draw_by_term("COX3", dat = glaTCdat, annotation = diamondl, species = "Gla")

draw_by_term("COX", dat = eveTCdat, annotation = diamondv, species = "Eve")


draw_by_term("glutathione S-transferase", dat = eveTCdat, annotation = diamondv, species = "Eve")


#DE GSTs
draw_by_id("TRINITY_DN376419_c0_g2_i3", dat = eveTCdat, annotation = diamondv, species = "Eve")
draw_by_id("TRINITY_DN376419_c0_g2_i6", dat = eveTCdat, annotation = diamondv, species = "Eve")

#E. cyaneus
draw_by_id("TRINITY_DN484746_c2_g1_i7", dat = ecyTCdat, annotation = diamondc, species = "Ecy")
draw_by_id("TRINITY_DN493696_c0_g1_i5", dat = ecyTCdat, annotation = diamondc, species = "Ecy")
draw_by_id("TRINITY_DN489637_c0_g2_i2", dat = ecyTCdat, annotation = diamondc, species = "Ecy")


draw_by_term("ribosomal protein", dat = eveTCdat, annotation = diamondv, species = "Eve")


draw_by_term("XP_015912981.1 heat shock protein 70 B2-like", 
             dat = eveALLdat, annotation = diamondv, species = "Eve")


draw_by_term("heat shock protein 70 B2-like", 
             dat = ecyALLdat, annotation = diamondc, species = "Ecy")

draw_by_term("XP_015912981.1 heat shock protein 70 B2-like", 
             dat = ecyALLdat, annotation = diamondc, species = "Gla")


draw_by_term("AEC33277.1 14-3-3 protein, partial", 
             dat = eveALLdat, annotation = diamondv, species = "Eve")
