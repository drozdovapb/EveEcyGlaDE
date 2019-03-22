#in case of transparency problems
#R should be compiled with the --use-cairo option
options(bitmapType = "cairo")

getTaxonomy <- function(thisassembly) {
    
  ##remember to choose
#  setwd("/media/main/sandbox/drozdovapb/DE/annotation/")
#  setwd("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/")
  setwd("/run/media/polina/Elements/transcriptome/DE/3-annotation")
    ##cut -f 13  GlaBCd6T24T_deep_FR.diamond.tsv >GlaBCd6T24T_deep_FR.taxid
    ## cut -f 13  GlaBCdTP1_cor.diamond.tsv >GlaBCdTP1_cor.taxid
    ##cut -f 13  EcyBCdTP1_cor.diamond.tsv >EcyBCdTP1_cor.taxid
    ##cut -f 13  GlaBCdTP1_cor2.diamond.tsv >GlaBCdTP1_cor2.taxid
    ids <- readLines(thisassembly)
    #remove multiple numbers in one entry, as they cannot be processed by NCBI
    sids <- sapply(strsplit(ids, split=";"), '[', 1)
    tids <- sort(table(sids))
    write.csv(tids, paste0(thisassembly, "tids.csv"))} 

matchTaxonomy <- function(tax_report, thisassembly) {
    #https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
    #got full taxonomy report
    #taxonomy report:
    trec <- read.delim(tax_report)
    tl <- strsplit(as.character(trec$lineage), " ")
    kingdom <- lapply(tl, function(x) x[length(x)-3])
    vkingdom <- unlist(lapply(kingdom, "[", 1))
    #qq <- cbind(as.numeric(vkingdom), as.numeric(tids))
    tids <- read.csv(paste0(thisassembly, "tids.csv"))
    qq <- data.frame(kingdom=as.numeric(vkingdom), num=as.numeric(tids$Freq), taxon=tids$sids)
    #summ <- aggregate(num~kingdom, qq, sum)
    #summ[order(summ$num),]
    return(qq)}


volcanoplot <- function(merged, thiscolor, thisalpha = 0.3, thesecondi) {
    sig <- ifelse(abs(merged$log2FoldChange) > 3 & merged$padj < 0.001, thiscolor, "black")
    p <- ggplot(merged, aes(log2FoldChange, -log10(padj))) +
        geom_point(col=sig, alpha=thisalpha) + ylim(0, 50) + xlim(-15, 15) +
      #geom_point(col=sig) + ylim(0, 50) + xlim(-15, 15) +
        ggtitle(paste0(thispecies, thesecondi[2], "_vs_", thesecondi[1]))
    return(p)
}

select.metazoa <- function(tbl, qq) {
    #tbl$taxon <- tbl$V13
    evem <- merge(x = tbl, y=qq, by="taxon", all.x=T, all.y=F)
    evem_metazoa <- evem[evem$kingdom==33208 & !is.na(evem$kingdom),]
    return(evem_metazoa)
}


deselect.contaminants <- function(tbl, qq) {
    #tbl$taxon <- tbl$V13
    evem <- merge(x = tbl, y=qq, by="taxon", all.x=T, all.y=F)
    evem_decont <- evem[evem$kingdom==33208 | is.na(evem$kingdom),] #NA also appears here, so not blasted should be retained...
    return(evem_decont)
}


perform.de <- function(dir, thispecies, theseconditions) {
    #Setup
    setwd(dir)
    #read samples table
    samples <- read.csv("./samples.csv")
    #create a vector of paths to files
    files <- file.path(dir, "gc", samples$salmon.folder, "quant.sf")
    names(files) <- paste0("sample", 1:6)
    #check whether everything is fine and all the files exist
    all(file.exists(files))
    
    
    #now, get colors right
    palette <- c("Blues", "Greens", "Oranges")
    ##colors <- c("blue4", "green4", "orange4")
    ##colors <- c("#0072B2", "#009E73", "#D55E00") #should be the color blind friendly panel but....
    colors <- c("#56B4E9", "#007656", "#D55E00") #I see it. Coblis also does. blue green orange.
    species <- c("Ecy", "Eve", "Gla")
    
    thispalette <- palette[which(species ==  thispecies)]
    thiscolor <- colors[which(species == thispecies)]
    
    ###Transcripts need to be associated with gene IDs for gene-level summarization.
    ###should I be fine with transcript-level summarization?
    
    #txi <- tximport(files, type = "salmon", txOut=TRUE)
    #it would never work at my laptop
    
    tocompare <- which(samples$species==thispecies & samples$condition %in% theseconditions)
    txi <- tximport(files[tocompare], type = "salmon", txOut=TRUE)
    #should be read in several minutes...
    
    #now let us construct a DESeq2 object
    library(DESeq2)
    
    sampleTable <- data.frame(condition=samples$condition[tocompare])
    #rownames(sampleTable) <- colnames(txi$counts) #it is NA
    rownames(sampleTable) <- samples$sample[tocompare]
    colnames(txi$counts) <- samples$sample[tocompare]
    
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    #80 to 20 Mb It's great but how would I now match back? It would be quite complicated
    
    #differential expression! 
    dds <- DESeq(dds)
    res <- results(dds)
    res$log10padj <- -log10(res$padj)
    
    ##mcols(res)$description
    
    resOrdered <- res[order(res$log2FoldChange),]
    
    write.csv(resOrdered, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_ordered.csv"))
    
    ##and some info for me
    message("number of DE genes for absolute log FC threshold = 1")
    temp <- resOrdered[complete.cases(resOrdered$padj), ]
    print(nrow(temp[abs(temp$log2FoldChange) > 1 & temp$padj < 0.001,]))
    message("number of DE genes for absolute log FC threshold = 3")
    print(nrow(temp[abs(temp$log2FoldChange) > 3 & temp$padj < 0.001,]))
    
    
    #Transformed data for some playing around
    vsd <- vst(dds, blind=FALSE)
    plotPCA(vsd, intgroup=c("condition"))
    
    svg(filename = paste0(thispecies, theseconditions[2],"_PCA",".svg"), 
        width=3.5, height=3.5)
    pl <- plotPCA(vsd, intgroup=c("condition")) + 
      ggtitle(paste0(thispecies, theseconditions[2])) +
      theme_bw() +  scale_color_manual(values = c("black", "red")) 
    print(pl)
    dev.off()
    
    sampleDists <- dist(t(assay(vsd)))
    ###clustering
    library("RColorBrewer")
    library("pheatmap") 
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(vsd)
    colnames(sampleDistMatrix) <- NULL
    
    colors <- colorRampPalette( rev(brewer.pal(9, thispalette)))(255)
    
    svg(filename = paste0(thispecies,theseconditions[2], "_clust",".svg"), 
        width=6, height=4)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors, fontsize = 16)
    #ok, it looks fine
    dev.off()
    #return(resOrdered)
    
}


filterTaxon <- function(thispecies, theseconditions) {
  
  #now, get colors right
  palette <- c("Blues", "Greens", "Oranges")
  ##colors <- c("blue4", "green4", "orange4")
  ##colors <- c("#0072B2", "#009E73", "#D55E00") #should be the color blind friendly panel but....
  ##colors <- c("#56B4E9", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thispalette <- palette[which(species ==  thispecies)]
  thiscolor <- colors[which(species == thispecies)]
  #and get alpha
  thisalpha <- ifelse(thispecies == "Ecy", 0.5, 0.3)
  
  resOrdered <- read.csv(paste0(thispecies,theseconditions[2],"_ordered.csv"), stringsAsFactors = F)
  
  ##Merge with 'annotation'

  #Annotation (best blast hit + also taxonomical information)
    diamond <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"),
#    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=F, stringsAsFactors = F)
  names(diamond)[1] <- "gene"
  diamond$V13 <- sapply(strsplit(diamond$V13, split=";"), '[', 1)
  #names(res)[1] <- "gene"
  #resOrdered$gene <- row.names(resOrdered)
  #now as I've read those, it should be different
  resOrdered$gene <- resOrdered$X
  
  premerged <- merge(y = diamond[,c(1,13,14)], x=as.data.frame(resOrdered), all.x=T, all.y=F , by="gene")
  
  #Annotation made by FunctionAnnotator  
  fa <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, "_AnnotationTable.txt"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=T, stringsAsFactors = F)
  names(fa)[1] <- "gene"
  
  merged <- merge(y = fa[,c(1:3,7:9,22)], x=premerged, all.x=T, all.y=F , by="gene")
  
  #rename
  names(merged)[which(names(merged) == "V13")] <- "taxon"
  names(merged)[which(names(merged) == "V14")] <- "best.nr.hit.diamond"
  
  mergedOrdered <- merged[order(merged$log2FoldChange),]
  #remove the genes with unidentifiable padj, or they will colonize the downstream things
  mergedOrdered <- mergedOrdered[complete.cases(mergedOrdered$padj),]
  write.csv(mergedOrdered, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_annot.csv"))
  
  de <- mergedOrdered[abs(mergedOrdered$log2FoldChange) > 3 & mergedOrdered$padj < 0.001,]
  write.csv(de, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_annot_de.csv"))
  
  p1 <- volcanoplot(merged, thiscolor, thesecondi = theseconditions)
   ggsave(paste0(thispecies,theseconditions[2],".png"), 
            p1, width = 3, device="png", height = 3)
  

  ####qq <- getTaxonomy(thisassembly)
  ####qq should be counted outside of this function.
  
  setwd(dir)
  
  #diamond <- read.delim("../../EveBCd6T24T_deep_FR.diamond.tsv", head=F, stringsAsFactors = F)
  
  #33208
  merged.metazoa <- select.metazoa(mergedOrdered, qq)
  merged.metazoa <- merged.metazoa[order(merged.metazoa$log2FoldChange),]
  #merged.metazoa <- merged.metazoa[complete.cases(merged.metazoa$padj),] #shouldn't be important now
  write.csv(merged.metazoa, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_metazoa.csv"))
  
  merged.metazoa.de <- merged.metazoa[abs(merged.metazoa$log2FoldChange) > 3 & merged.metazoa$padj < 0.001,]    
  #merged.metazoa.de <- merged.metazoa.de[complete.cases(merged.metazoa.de),] #again, sholdn't make any diff now
  write.csv(merged.metazoa.de,
            paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_metazoa_de.csv"))
  
  p2 <- volcanoplot(merged.metazoa, thiscolor, thesecondi = theseconditions) + 
    ggtitle(paste0(thispecies, theseconditions[2], "_vs_", theseconditions[1], "_Metazoa"))
  #p
  ggsave(paste0(thispecies,theseconditions[2],"_metazoa.png"), 
         p2, device = "png", width = 3, height = 3)
  
  
  merged.decont <- deselect.contaminants(merged, qq)
  merged.decont <- merged.decont[order(merged.decont$log2FoldChange),]
  write.csv(merged.decont, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_decont.csv"))
  
  merged.decont.de <- merged.decont[abs(merged.decont$log2FoldChange) > 3 & merged.decont$padj < 0.001,]    
  merged.decont.de <- merged.decont.de[complete.cases(merged.decont.de$gene),] #but this is important!
  write.csv(merged.decont.de,
            paste0(thispecies,theseconditions[2],"_decont_de.csv"))
  
  p3 <- volcanoplot(merged.decont, thiscolor, thesecondi = theseconditions) + 
    ggtitle(paste0(thispecies, theseconditions[2], "_vs_", theseconditions[1], "-contam-n"))
  #p
  ggsave(paste0(thispecies,theseconditions[2],"_decont.png"), 
         p3, device = "png", width = 3, height = 3)
  
}


##Another function for filterTaxon
filterTaxonReduced <- function(thispecies, theseconditions, logFCthreshold = 3) {
  
  #now, get colors right
  palette <- c("Blues", "Greens", "Oranges")
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thispalette <- palette[which(species ==  thispecies)]
  thiscolor <- colors[which(species == thispecies)]
  #and get alpha
  thisalpha <- ifelse(thispecies == "Ecy", 0.5, 0.3)
  
  resOrdered <- read.csv(paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_ordered.csv"), stringsAsFactors = F)
  
  ##Merge with 'annotation'
  
  #Annotation (best blast hit + also taxonomical information)
  diamond <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    #paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"),
    paste0("~/Documents/transcriptome_annotation/", thisassembly, ".diamond.tsv"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=F, stringsAsFactors = F)
  names(diamond)[1] <- "gene"
  diamond$V13 <- sapply(strsplit(diamond$V13, split=";"), '[', 1)
  #names(res)[1] <- "gene"
  #resOrdered$gene <- row.names(resOrdered)
  #now as I've read those, it should be different
  resOrdered$gene <- resOrdered$X
  
  premerged <- merge(y = diamond[,c(1,13,14)], x=as.data.frame(resOrdered), all.x=T, all.y=F , by="gene")
  
  #Annotation made by FunctionAnnotator  
  fa <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    paste0("~/Documents/transcriptome_annotation/", thisassembly, "_AnnotationTable.txt"),
    #paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, "_AnnotationTable.txt"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=T, stringsAsFactors = F)
  names(fa)[1] <- "gene"
  
  merged <- merge(y = fa[,c(1:3,7:9,22)], x=premerged, all.x=T, all.y=F , by="gene")
  
  #rename
  names(merged)[which(names(merged) == "V13")] <- "taxon"
  names(merged)[which(names(merged) == "V14")] <- "best.nr.hit.diamond"
  
  mergedOrdered <- merged[order(merged$log2FoldChange),]
  #remove the genes with unidentifiable padj, or they will colonize the downstream things
  mergedOrdered <- mergedOrdered[complete.cases(mergedOrdered$padj),]
  write.csv(mergedOrdered, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1], "_annot.csv"))
  
  de <- mergedOrdered[abs(mergedOrdered$log2FoldChange) > logFCthreshold & mergedOrdered$padj < 0.001,]
  write.csv(de, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1], logFCthreshold, "_annot_de.csv"))
  
  p1 <- volcanoplot(merged, thiscolor, thesecondi = theseconditions)
  ggsave(paste0(thispecies,theseconditions[2], theseconditions[1], logFCthreshold,".png"), 
         p1, width = 3, device="png", height = 3)
  
}



##Another function for filterTaxon
filterTaxonReduced <- function(thispecies, theseconditions, logFCthreshold = 3) {
  
  #now, get colors right
  palette <- c("Blues", "Greens", "Oranges")
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thispalette <- palette[which(species ==  thispecies)]
  thiscolor <- colors[which(species == thispecies)]
  #and get alpha
  thisalpha <- ifelse(thispecies == "Ecy", 0.5, 0.3)
  
  resOrdered <- read.csv(paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1],"_ordered.csv"), stringsAsFactors = F)
  
  ##Merge with 'annotation'
  
  #Annotation (best blast hit + also taxonomical information)
  diamond <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=F, stringsAsFactors = F)
  names(diamond)[1] <- "gene"
  diamond$V13 <- sapply(strsplit(diamond$V13, split=";"), '[', 1)
  #names(res)[1] <- "gene"
  #resOrdered$gene <- row.names(resOrdered)
  #now as I've read those, it should be different
  resOrdered$gene <- resOrdered$X
  
  premerged <- merge(y = diamond[,c(1,13,14)], x=as.data.frame(resOrdered), all.x=T, all.y=F , by="gene")
  
  #Annotation made by FunctionAnnotator  
  fa <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, "_AnnotationTable.txt"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=T, stringsAsFactors = F)
  names(fa)[1] <- "gene"
  
  merged <- merge(y = fa[,c(1:3,7:9,22)], x=premerged, all.x=T, all.y=F , by="gene")
  
  #rename
  names(merged)[which(names(merged) == "V13")] <- "taxon"
  names(merged)[which(names(merged) == "V14")] <- "best.nr.hit.diamond"
  
  mergedOrdered <- merged[order(merged$log2FoldChange),]
  #remove the genes with unidentifiable padj, or they will colonize the downstream things
  mergedOrdered <- mergedOrdered[complete.cases(mergedOrdered$padj),]
  write.csv(mergedOrdered, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1], "_annot.csv"))
  
  de <- mergedOrdered[abs(mergedOrdered$log2FoldChange) > logFCthreshold & mergedOrdered$padj < 0.001,]
  write.csv(de, paste0(thispecies,theseconditions[2],"_vs_",theseconditions[1], logFCthreshold, "_annot_de.csv"))
  
  p1 <- volcanoplot(merged, thiscolor, thesecondi = theseconditions)
  ggsave(paste0(thispecies,theseconditions[2], theseconditions[1], logFCthreshold,".png"), 
         p1, width = 3, device="png", height = 3)
  
}


##interspecific
perform.de.inter <- function(dir, thispecies, theseconditions) {
  #Setup
  setwd(dir)
  #read samples table
  samples <- read.csv("./samples.csv")
  #create a vector of paths to files
  files <- file.path(dir, "gc", samples$salmon.folder, "quant.sf")
  names(files) <- paste0("sample", 1:6)
  #check whether everything is fine and all the files exist
  all(file.exists(files))
  
  
  #now, get colors right
  #palette <- c("Blues", "Greens", "Oranges")
  ##colors <- c("blue4", "green4", "orange4")
  ##colors <- c("#0072B2", "#009E73", "#D55E00") #should be the color blind friendly panel but....
  colors <- c("#56B4E9", "#007656", "#D55E00") #I see it. Coblis also does. blue green orange.
  species <- c("Ecy", "Eve", "Gla")
  
  #thispalette <- palette[which(species ==  thispecies)]
  #thiscolor <- colors[which(species == thispecies)]
  thiscolor <- "#FF0000"
  
  ###Transcripts need to be associated with gene IDs for gene-level summarization.
  ###should I be fine with transcript-level summarization?
  
  #txi <- tximport(files, type = "salmon", txOut=TRUE)
  #it would never work at my laptop
  
  tocompare <- which(samples$species %in% thispecies & samples$condition_ %in% theseconditions)
  txi <- tximport(files[tocompare], type = "salmon", txOut=TRUE)
  #should be read in several minutes...
  
  #now let us construct a DESeq2 object
  library(DESeq2)
  
  sampleTable <- data.frame(condition=samples$condition_[tocompare], species=samples$species[tocompare])
  #rownames(sampleTable) <- colnames(txi$counts) #it is NA
  rownames(sampleTable) <- samples$sample[tocompare]
  colnames(txi$counts) <- samples$sample[tocompare]
  
  dds <- DESeqDataSetFromTximport(txi, sampleTable, ~species)
  
  ## this time don't get rid of anything, we have way too less data...
  #keep <- rowSums(counts(dds)) >= 10
  #dds <- dds[keep,]
  #80 to 20 Mb It's great but how would I now match back? It would be quite complicated
  
  #differential expression! 
  dds <- DESeq(dds)
  res <- results(dds)
  res$log10padj <- -log10(res$padj)
  
  ##mcols(res)$description
  
  #Transformed data for some playing around
  #vsd <- vst(dds, blind=FALSE)
  #plotPCA(vsd, intgroup=c("species"))
  
  #svg(filename = paste0(thispecies[1], thispecies[2],"_PCA",".svg"), 
  #    width=3.5, height=3.5)
  #pl <- plotPCA(vsd, intgroup=c("species")) + 
  #  ggtitle(paste0(thispecies, theseconditions[2])) +
  #  theme_bw() +  scale_color_manual(values = c("black", "red")) 
  #print(pl)
  #dev.off()
  
  # sampleDists <- dist(t(assay(vsd)))
  # ###clustering
  # library("RColorBrewer")
  # library("pheatmap") 
  # sampleDistMatrix <- as.matrix(sampleDists)
  # rownames(sampleDistMatrix) <- colnames(vsd)
  # colnames(sampleDistMatrix) <- NULL
  # 
  # colors <- colorRampPalette( rev(brewer.pal(9, thispalette)))(255)
  # 
  # # svg(filename = paste0(thispecies,theseconditions[2], "_clust",".svg"), 
  # #     width=6, height=4)
  #  pheatmap(sampleDistMatrix,
  #           clustering_distance_rows=sampleDists,
  #           clustering_distance_cols=sampleDists,
  #           col=colors, fontsize = 16)
  # # #ok, it looks fine
  # # dev.off()
  # # 
  resOrdered <- res[order(res$log2FoldChange),]
  
  write.csv(resOrdered, paste0(thispecies[1], thispecies[2], theseconditions[1], "_ordered.csv"))
  
  
  return(resOrdered)
  
}


filterTaxonInter <- function(thispecies, theseconditions, logFCthreshold = 3) {
  
  #now, get colors right
  palette <- c("Blues", "Greens", "Oranges")
  ##colors <- c("blue4", "green4", "orange4")
  ##colors <- c("#0072B2", "#009E73", "#D55E00") #should be the color blind friendly panel but....
  ##colors <- c("#56B4E9", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  colors <- c("#4087AF", "#007656", "#D55E00") #I see it. Coblis does as well. Blue, green, orange
  species <- c("Ecy", "Eve", "Gla")
  
  thispalette <- palette[which(species ==  thispecies[2])]
  thiscolor <- colors[which(species == thispecies[2])]
  #and get alpha
  thisalpha <- ifelse(thispecies == "Ecy", 0.5, 0.3)
  
  #resOrdered <- read.csv(paste0(thispecies[1],thispecies[2],theseconditions,"_ordered.csv"), stringsAsFactors = F)
  ##or 
  resOrdered <- read.csv(paste0(thispecies[1],thispecies[2],theseconditions[2],"_ordered.csv"), stringsAsFactors = F)
  
  #Annotation (best blast hit + also taxonomical information)
  diamond <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    #paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"),
    paste0("~/Documents/transcriptome_annotation/", thisassembly, ".diamond.tsv"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=F, stringsAsFactors = F)
  names(diamond)[1] <- "gene"
  diamond$V13 <- sapply(strsplit(diamond$V13, split=";"), '[', 1)
  #names(res)[1] <- "gene"
  #resOrdered$gene <- row.names(resOrdered)
  #now as I've read those, it should be different
  resOrdered$gene <- resOrdered$X
  
  premerged <- merge(y = diamond[,c(1,13,14)], x=as.data.frame(resOrdered), all.x=T, all.y=F , by="gene")
  
  #Annotation made by FunctionAnnotator  
  fa <- read.delim(
    #        paste0("/media/main/sandbox/drozdovapb/DE/annotation/", thisassembly, ".diamond.tsv"), 
    #paste0("/run/media/polina/Elements/transcriptome/DE/3-annotation/", thisassembly, "_AnnotationTable.txt"),
    paste0("~/Documents/transcriptome_annotation/", thisassembly, "_AnnotationTable.txt"),
    #    paste0("/media/drozdovapb/Elements/transcriptome/DE/3-annotation/", thisassembly, ".diamond.tsv"), 
    head=T, stringsAsFactors = F)
  names(fa)[1] <- "gene"
  
  merged <- merge(y = fa[,c(1:3,7:9,22)], x=premerged, all.x=T, all.y=F , by="gene")
  
  #rename
  names(merged)[which(names(merged) == "V13")] <- "taxon"
  names(merged)[which(names(merged) == "V14")] <- "best.nr.hit.diamond"
  
  mergedOrdered <- merged[order(merged$log2FoldChange),]
  ##remove the genes with unidentifiable padj, or they will colonize the downstream things
  mergedOrdered <- mergedOrdered[complete.cases(mergedOrdered$padj),]
  write.csv(mergedOrdered, paste0(thispecies[1],thispecies[2],theseconditions[2],"_vs_",theseconditions[1], "_annot.csv"))
  
  de <- mergedOrdered[abs(mergedOrdered$log2FoldChange) > logFCthreshold & mergedOrdered$padj < 0.001,]
  write.csv(de, paste0(thispecies[1],thispecies[2],theseconditions[2],"_vs_",theseconditions[1], 
                       logFCthreshold, "_annot_de.csv"))
  
  p1 <- volcanoplot(merged, thiscolor, thesecondi = theseconditions)
  ggsave(paste0(thispecies[1],thispecies[2],theseconditions[2], theseconditions[1], logFCthreshold,".png"), 
         p1, width = 3, device="png", height = 3)
  
}
