setwd("/media/DATA/учебное-лабораторное/лаборатория/R и прочее/термостабильность белков")


evd <- read.table("./7-protein-lengths/EveBCdT_deep_decontaminated.fasta.transdecoder.pep.sl.complete.fasta.size", stringsAsFactors=FALSE) # DEEP
eve <- read.table("./7-protein-lengths/EveBCdTP1_cor_decontaminated.fasta.transdecoder.pep.sl.complete.fasta.size", stringsAsFactors=FALSE)
ecy <- read.table("./7-protein-lengths/EcyBCdTP1_cor_decontaminated.fasta.transdecoder.pep.sl.complete.fasta.size", stringsAsFactors=FALSE)
gla <- read.table("./7-protein-lengths/GlaBCdTP1_cor_decontaminated.fasta.transdecoder.pep.sl.complete.fasta.size", stringsAsFactors=FALSE)

boxplot(list(evd[,2], eve[,2], ecy[,2], gla[,2]), range=20)

median(evd[,2])
median(eve[,2])
median(ecy[,2])
median(gla[,2])

mean(evd[,2])
mean(eve[,2])
mean(ecy[,2])
mean(gla[,2])

wilcox.test(evd[,2], eve[,2]) # две сборки сильно отличаются
wilcox.test(eve[,2], ecy[,2]) # байкальские одинаковы
wilcox.test(eve[,2], gla[,2]) # лакустрис отличается от
wilcox.test(ecy[,2], gla[,2]) # обоих байкальских





#Cellular Proteomes Have Broad Distributions of Protein Stability, 2010
#calculate delta G
dg <- function(l, t) {
  s <- (-5.03*l - 41.6) + (-0.062*l + 0.53)*(t - 373.5) - t*(-0.0168*l - 0.085) - t*(-0.062*l + 0.53)*log(t/385)
  return(s)
}

gd6 <- sapply(evd[,2], dg, 280) #Eve deep
ge6 <- sapply(eve[,2], dg, 280) #Eve
gc6 <- sapply(ecy[,2], dg, 280) #Ecy
gg6 <- sapply(gla[,2], dg, 280) #Gla

gd25 <- sapply(evd[,2], dg, 298)
ge25 <- sapply(eve[,2], dg, 298)
gc25 <- sapply(ecy[,2], dg, 298)
gg25 <- sapply(gla[,2], dg, 298)

gd37 <- sapply(evd[,2], dg, 310)
ge37 <- sapply(eve[,2], dg, 310)
gc37 <- sapply(ecy[,2], dg, 310)
gg37 <- sapply(gla[,2], dg, 310)

boxplot(list(gd6, gd25), range=20)

g6 <- list(gd6, ge6, gc6, gg6)
g25 <- list(gd25, ge25, gc25, gg25)

sapply(g6, median)
sapply(g25, median)





# fraction of proteins of length l that are unfolded at temperature t
nu <- function(g, t) {
  nu <- 1/(1 + exp(-g/(0.083*t)))
  return(nu)
}

nd6 <- sapply(gd6, nu, 280)
ne6 <- sapply(ge6, nu, 280)
nc6 <- sapply(gc6, nu, 280)
ng6 <- sapply(gg6, nu, 280)
c(sum(nd6 > 0.3), sum(ne6 > 0.3), sum(nc6 > 0.3), sum(ng6 > 0.3))

nd25 <- sapply(gd25, nu, 298)
ne25 <- sapply(ge25, nu, 298)
nc25 <- sapply(gc25, nu, 298)
ng25 <- sapply(gg25, nu, 298)
c(sum(nd25 > 0.3), sum(ne25 > 0.3), sum(nc25 > 0.3), sum(ng25 > 0.3))

median(nd25)
median(ne25)
median(nc25)
median(ng25)

nd37 <- sapply(gd37, nu, 310)
ne37 <- sapply(ge37, nu, 310)
nc37 <- sapply(gc37, nu, 310)
ng37 <- sapply(gg37, nu, 310)
c(sum(nd37 > 0.3), sum(ne37 > 0.3), sum(nc37 > 0.3), sum(ng37 > 0.3))
c(length(nd37), length(ne37), length(nc37), length(ng37))
c(sum(nd37 > 0.3)/length(nd37), sum(ne37 > 0.3)/length(ne37), sum(nc37 > 0.3)/length(nc37), sum(ng37 > 0.3)/length(ng37))


