dir <- "/media/drozdovapb/big/Research/Projects/DE/"; setwd(dir)

source("salmon/de_functions_new+decontamination.R")


# get_decon_list <- function(thisassembly) {
#     getTaxonomy(paste0(thisassembly, ".taxid"))
#     qq <- matchTaxonomy(paste0(thisassembly, "_taxreport.txt"),
#                         paste0(thisassembly, ".taxid"))
#     
#     #get back to the target working directory
#     #dir <- "/media/drozdovapb/big/Research/Projects/DE/salmon/eve_deep"; setwd(dir)
#     diamond <- read.delim(
#         paste0("/media/drozdovapb/big/Research/Projects/DE/annotation/", 
#                thisassembly, ".diamond.tsv"), 
#         head=F, stringsAsFactors = F)
#     
#     diamond$taxon <- diamond$V13
#     wtaxa <- merge(x = diamond, y=qq, by="taxon", all.x=T, all.y=F)
#     wtaxa_decont <- wtaxa[wtaxa$kingdom==33208 | is.na(wtaxa$kingdom),] #NA also appears here, so not blasted should be retained...
#     writeLines(wtaxa_decont$V1, paste0(thisassembly, "_names_decontaminated.txt"))
#     #and also some workaround for sequences that were not recognized by diamond!!
# }

#We'll go the other way round
get_contaminant_list <- function(thisassembly) {
    getTaxonomy(paste0(thisassembly, ".taxid"))
    qq <- matchTaxonomy(paste0(thisassembly, "_taxreport.txt"),
                        paste0(thisassembly, ".taxid"))
    #get back to the target working directory
    #dir <- "/media/drozdovapb/big/Research/Projects/DE/salmon/eve_deep"; setwd(dir)
    diamond <- read.delim(
        paste0("/media/drozdovapb/big/Research/Projects/DE/annotation/", 
               thisassembly, ".diamond.tsv"), 
        head=F, stringsAsFactors = F)
    
    diamond$taxon <- as.numeric(diamond$V13)
    wtaxa <- merge(x = diamond, y=qq, by="taxon", all.x=T, all.y=F)
    #wtaxa_cont1 <- wtaxa[!(wtaxa$kingdom == 33208) & !is.na(wtaxa$kingdom),]
    #wtaxa_metazoa <- wtaxa[(wtaxa$kingdom == 33208) & !is.na(wtaxa$kingdom),]
    #I have a great problem with NAs!
    wtaxa[is.na(wtaxa$kingdom),]$kingdom <- "unknown"
    wtaxa_cont <- wtaxa[!(wtaxa$kingdom == 33208),]
    wtaxa_cont$V1 <- paste0(">", wtaxa_cont$V1)
    writeLines(wtaxa_cont$V1, paste0(thisassembly, "_contamination.txt"))
    #and also some workaround for sequences that were not recognized by diamond!!
    return(wtaxa_cont$V1)
}

setwd(dir)

decontaminate <- function(thisassembly, assembly_cont) {
    setwd(dir)
    assembly_all <- readLines(paste0(thisassembly, ".names"))
    assembly_all_short <- sapply(" .*", gsub, "", assembly_all)
    assembly_decontaminated <- setdiff(assembly_all_short, assembly_cont)
    assembly_decontaminated <- sapply(">", gsub, "", assembly_decontaminated)
    writeLines(assembly_decontaminated, paste0("./", thisassembly, "_decont.txt"))
}



#thisassembly <- "EcyBCdTP1_cor"
ecy_cont <- get_contaminant_list("EcyBCdTP1_cor")
#146255... does it sound fine? Yes, it does
#grep \> EcyBCdTP1_cor.fasta >EcyBCdTP1_cor.names
decontaminate("EcyBCdTP1_cor", ecy_cont)


eve_cont <- get_contaminant_list("EveBCdTP1_cor")
#grep \> EveBCdTP1_cor.fasta >EveBCdTP1_cor.names
#grep \> GlaBCdTP1_cor.fasta >GlaBCdTP1_cor.names
decontaminate("EveBCdTP1_cor", eve_cont)



#finally!

thisassembly <- "GlaBCdTP1_cor"
gla_cont <- get_contaminant_list("GlaBCdTP1_cor")
decontaminate("GlaBCdTP1_cor", gla_cont)

thisassembly <- "EveBCd6T24T_deep_FR"
eve_deep_cont <- get_contaminant_list("EveBCd6T24T_deep_FR")

#stopped_here!!!
#grep \> annotation/EveBCd6T24T_deep_FR.fasta >|EveBCd6T24T_deep_FR.names
decontaminate("EveBCd6T24T_deep_FR", eve_deep_cont)
#eve_deep_cont <- get_contaminant_list("EveBCd6T24T_deep_FR")

