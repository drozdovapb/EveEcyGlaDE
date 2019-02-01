#[polina@ale transcriptome_annotation]$ getorf EcyCd24up_unknown.fasta -outseq EcyCd24up_unknown.faa


#install.packages("seqinr")
library(seqinr)

setwd("/homes/brauerei/polina/Documents/transcriptome_annotation")

seqs <- read.fasta("./EcyCd24up_unknown.faa")
#seqs <- read.fasta("~/Documents/Paper1_stresses/crust_MTs.fasta")

#seq1 <- seqs[[1]]

##get length
#length(seq1)
##get aa distribution
#ts <- table(seq1)             
##which are aromatic? WYF? H?
#aromatic <- which(names(ts)=="w" | names(ts) == "y" | names(ts) == "f" | names(ts) == "h")
#percent_aromatic <- sum(ts[aromatic])/length(seq1) # eg 15% in here
#cystein <- which(names(ts)=="c")
#percent_cys <- ts[cystein]/length(seq1)

#with apply? 

#starts with methionine
#seqs <- seqs[sapply(seqs, FUN = function(x) x[[1]] == "m")]

lengths <- sapply(seqs, length)

common_table <- data.frame(lengths = lengths, aromatic = NA, cystein = NA)

for (i in 1:length(seqs)) {
  ts <- table(seqs[[i]])             
  #which are aromatic? WYF? H?
  aromatic <- which(names(ts)=="w" | names(ts) == "y" | names(ts) == "f" | names(ts) == "h")
  percent_aromatic <- sum(ts[aromatic])/length(seqs[[i]]) * 100 # eg 15% in here
  common_table$aromatic[i] <- percent_aromatic
  #cystein is very importna!
  cystein <- which(names(ts)=="c")
  if (!length(cystein)) percent_cys <- 0
  if (length(cystein)) percent_cys <- ts[cystein]/length(seqs[[i]]) * 100
  common_table$cystein[i] <- percent_cys
}


cfilt <- common_table[common_table$aromatic < 10 & common_table$cystein > 10 & common_table$lengths > 25 & 
                        common_table$lengths < 100,]
#passed by most known (64/66) crustacean MT proteins...
#then let's assume it's a good criterion

#158 predicted protein sequences
#how many contigs? 
unique(substr(row.names(cfilt), 1, 25)) #all of them are unique
# 50 unique contigs

# well... what else can I do? 
# starting with M? 

cfilt



setwd("/homes/brauerei/polina/Documents/Paper1_stresses/MTs")
v <- readLines("EveCd24_upregs.txt")
l <- readLines("EcyCd24_upregs.txt")
ivl <- intersect(v, l)
#already just a little
#wait, am I using the same assembly?
intersect(row.names(cfilt), ivl) #nothing


intersect(substr(row.names(common_table), 1, 25), ivl) #nothing

intersect(substr(row.names(cfilt), 1, 25), ivl) #nothing



#nope. Maybe look for patterns?
#and relax the conditions?

c <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/ecy/EcyCd03_ordered.csv")
c24 <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/ecy/EcyCd24_ordered.csv")
vc <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/ecy/EveCd03_ordered.csv")
gc <- read.csv("/run/media/polina/Elements/transcriptome/DE/4-salmon/ecy/GlaCd03_ordered.csv")
  
c <- c[c$log10padj < 0.05 & c$log2FoldChange > 1, ]
vc <- vc[vc$log10padj < 0.05 & vc$log2FoldChange > 1, ]
lookfor <- intersect(c$X, vc$X)
gc <- gc[gc$log10padj < 0.05 & gc$log2FoldChange > 1, ]
c24 <- c24[c24$log10padj < 0.05 & c24$log2FoldChange > 1,]


i <- intersect(vc$X, intersect(c$X, c24$X))
writeLines(i[-1], "EveEcy_ids.txt")

#[polina@ale MTs]$ samtools faidx ../../transcriptome_annotation/EcyBCdTP1_cor.fasta < EveEcy_ids.txt > Ecy_putmts.fa


library(seqinr)

seqs <- read.fasta("trinity_dn476957_c0_g1_i2.orf", as.string = T, forceDNAtolower = F)


seqs <- read.fasta("trinity_dn476957_c0_g1_i2.orf")
seqs <- seqs[sapply(seqs, FUN = function(x) x[[1]] == "m")]

lengths <- sapply(seqs, length)


common_table <- data.frame(lengths = lengths, aromatic = NA, cystein = NA)


for (i in 1:length(seqs)) {
  ts <- table(seqs[[i]])             
  #which are aromatic? WYF? H?
  aromatic <- which(names(ts)=="w" | names(ts) == "y" | names(ts) == "f" | names(ts) == "h")
  percent_aromatic <- sum(ts[aromatic])/length(seqs[[i]]) * 100 # eg 15% in here
  common_table$aromatic[i] <- percent_aromatic
  #cystein is very importna!
  cystein <- which(names(ts)=="c")
  if (!length(cystein)) percent_cys <- 0
  if (length(cystein)) percent_cys <- ts[cystein]/length(seqs[[i]]) * 100
  common_table$cystein[i] <- percent_cys
}

common_table[common_table$cystein > 10 & common_table$aromatic < 10, ]
