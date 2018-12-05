## This script produces contig names attributed to a particular taxon

##set wd
wdir <- "/homes/brauerei/polina/Documents/transcriptome_annotation/taxonomy/"
##also make sure there is the ../names folder
setwd(wdir)


#######################FUNCTIONS###################################################################

##this function summarizes and outputs unique IDs for each assembly.
##Actually, it doesn't need to be done for each assembly but I'm lazy to make a comprehensive list.
getTaxonomy <- function(thisassembly) {
  ids <- readLines(paste0(thisassembly, ".taxid"))
  #remove multiple numbers in one entry, as they cannot be processed by NCBI
  sids <- sapply(strsplit(ids, split=";"), '[', 1)
  tids <- sort(table(sids))
  write.csv(tids, paste0(thisassembly, "_tids.csv"))
  ##These _tids.csv files contain unique ids + their frequence
} 

## This function reads taxonomy reports from NCBI and outputs a table of report -- taxon
matchTaxonomy <- function(thisassembly) {
  ## Read taxonomy report file
  trec <- read.delim(paste0(thisassembly, "_taxreport.txt"))
  tl <- strsplit(as.character(trec$lineage), " ")
  ## 3 is for phylum level. Choose whatever suits best
  kingdom <- lapply(tl, function(x) x[length(x)-3])
  vkingdom <- unlist(lapply(kingdom, "[", 1))
  tids <- read.csv(paste0(thisassembly, "_tids.csv"))
  qq <- data.frame(kingdom=as.numeric(vkingdom), num=as.numeric(tids$Freq), taxon=tids$sids)
  ## so, this is the conversion table we need
  write.csv(qq, paste0(thisassembly, "_taxid_conversion.csv"))
  ## in addition, we can report & write most abundant taxa (top 5)
  summ <- aggregate(num~kingdom, qq, sum)
  message("33208 = Metazoa \n 5878 == Ciliophora")
  head(summ[order(summ$num, decreasing = T),], 5)
  }

## This function selects contigs with various annotations
select.everything <- function(thisassembly) {
  ## read diamond annotations
  tbl <- read.delim(paste0("../", thisassembly, ".diamond.tsv"), head = F, stringsAsFactors = F)
  tbl$taxon <- tbl$V13
  ## read conversion table
  qq <- read.csv(paste0(thisassembly, "_taxid_conversion.csv"))
  evem <- merge(x = tbl, y=qq, by="taxon", all.x=T, all.y=F)
  ## Select Metazoa
  evem_metazoa <- evem[evem$kingdom==33208 & !is.na(evem$kingdom), "V1"]
  writeLines(evem_metazoa, paste0("../names/", thisassembly, "_Metazoa.names.txt"))
  ## Select Ciliophora
  evem_ciliophora <- evem[evem$kingdom==5878 & !is.na(evem$kingdom), "V1"]
  writeLines(evem_ciliophora, paste0("../names/", thisassembly, "_Ciliophora.names.txt"))
  ## Select other stuff
  ## Here we use another logic (all remaining), or otherwise we get weird results
  evem_other <- setdiff(evem$V1, evem[evem$kingdom==5878 | evem$kingdom==33208, "V1",])
  writeLines(evem_other, paste0("../names/", thisassembly, "_Other.names.txt"))
}


write.names <- function(assembly) {
  ##read diamond file
  tbl <- read.delim(paste0(wdir, thisassembly, ".diamond.tsv"), head=F, stringsAsFactors = F)
}

####################INSTRUCTIONS FOR THE MAIN PART#################################################

##First step, faster with bash
##cd taxonomy ##mkdir taxonomy first if it doesn't exist
##cut -f 13  EveBCdT_deep.diamond.tsv >EveBCdT_deep.taxid
##cut -f 13  EveBCdTP1_cor.diamond.tsv >EveBCdTP1_cor.taxid
##cut -f 13  EcyBCdTP1_cor.diamond.tsv >EcyBCdTP1_cor.taxid
##cut -f 13  GlaBCdTP1_cor2.diamond.tsv >GlaBCdTP1_cor2.taxid

##Our next problem is that NCBI taxonomy portal wouldn't handle 200k numbers
##For this we use getTaxonomy, which summarizes 

getTaxonomy("EveBCdT_deep")
getTaxonomy("EveBCdTP1_cor")
getTaxonomy("EcyBCdTP1_cor")
getTaxonomy("GlaBCdTP1_cor2")

##These _tids.csv files contain unique ids + their frequence
##Only IDs need to be copied into the form
##https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
##These files need to be saved as <thisassembly>_taxreport.txt
##Yes, I did it manually. 
##It's easier to save 4 files than to write up some code for API.

matchTaxonomy("EveBCdT_deep")
matchTaxonomy("EveBCdTP1_cor")
matchTaxonomy("EcyBCdTP1_cor")
matchTaxonomy("GlaBCdTP1_cor2")

##Now we're ready to select contigs
select.everything("EveBCdT_deep")
select.everything("EveBCdTP1_cor")
select.everything("EcyBCdTP1_cor")
select.everything("GlaBCdTP1_cor2")
