options(stringsAsFactors = F)
library(UpSetR)
library(ggplot2)

## re = rsem /edger;
## rd = rsem / deseq2; 
## se = salmon / edgr;
## sd = salmon / deseq2

## later make it a function
comparison <- "PB24_vs_B24h"
species <- "Ecy"

## a small function to extract DE genes
threshold_p <- 0.05
selectDEd <- function(tbl) {
  tbl[(tbl$log2FoldChange > 1 | tbl$log2FoldChange < -1) & tbl$padj < threshold_p, ]
}
selectDEe <- function(tbl) {
  tbl[(tbl$logFC > 1 | tbl$logFC < -1) & tbl$FDR < threshold_p, ]
}

get_common_list <- function(comparison, species) {
      ## I have 4 files to compare, just need to read them all and combine
    est_method <- c("salmon", "salmon", "bowtie_rsem", "bowtie_rsem")
    de_method <- c("DESeq2", "edgeR", "DESeq2", "edgeR")
    ## paths to 4 files
    paths <- paste0(est_method, "/", species, ".isoform.counts.matrix.", 
                    comparison, ".", de_method, ".DE_results")
    ## read files
    ## salmon / DESeq2
    sd <- read.delim(paths[1])
    ## salmon / edgeR
    se <- read.delim(paths[2])
    ## bowtie_rsem / DESeq2
    rd <- read.delim(paths[3])
    ## bowtie_rsem / edgeR
    re <- read.delim(paths[4])
  
    #select DE genes...
    sd_DE <- row.names(selectDEd(sd))
    se_DE <- row.names(selectDEe(se))
    rd_DE <- row.names(selectDEd(rd))
    re_DE <- row.names(selectDEe(re))
  
    ## intersections at least
    ## upSet figures
    lI <- list(salmon.DESeq2 = sd_DE, salmon.edgeR = se_DE, 
               rsem.DESeq2 = rd_DE, rsem.edgeR = re_DE)
    upset(fromList(lI), nintersects = NA)
    ## get the list of common elements
    common_list <- Reduce(intersect, lI)
    
    sd_common <- sd[row.names(sd) %in% common_list, ]
    se_common <- se[row.names(se) %in% common_list, ]
    rd_common <- rd[row.names(rd) %in% common_list, ]
    re_common <- re[row.names(re) %in% common_list, ]
    
    result <- data.frame(Gene = common_list,
                         salmon.DESeq2.FC = sd_common$log2FoldChange, 
                         salmon.edgeR.FC = se_common$logFC,
                         rsem.DESeq2.FC = rd_common$log2FoldChange, 
                         rsem.edgeR.FC = re_common$logFC)
## if it works, then add p-values and so on!
    # the line below: write just once
    common_DE <- diamond[diamond$V1 %in% common_list, ]
    write.csv(common_DE, paste0(comparison, ".", species, "common_DE.csv"))
    return(result)
}

## todo: this depends on the species. Think about it
diamond <- read.delim("~/Research/Projects/DE/annotation/EcyBCdTP1_cor.diamond.tsv", head = F)
ecy_a03 <- get_common_list("PB03_vs_B03", "Ecy")
ecy_a24 <- get_common_list("PB24_vs_B24h", "Ecy")
ecy_p03 <- get_common_list("Ph03_vs_PB03", "Ecy")
ecy_p24 <- get_common_list("Ph24_vs_PB24", "Ecy")
ecy_ap03 <- get_common_list("Ph03_vs_B03", "Ecy")
ecy_ap24 <- get_common_list("Ph24_vs_B24h", "Ecy")



## I'm afraid but need to do it. What's the correlation?


  mergedDE <- merge(ecy_a24, ecy_p24, by = "Gene", all = T)
  
  xmin <- min(mergedDE$salmon.DESeq2.FC.x, na.rm = T)
  xmax <- max(mergedDE$salmon.DESeq2.FC.x,  na.rm = T)
  ymin <- min(mergedDE$salmon.DESeq2.FC.y,  na.rm = T)
  ymax <- max(mergedDE$salmon.DESeq2.FC.y,  na.rm = T)
  fmin <- min(xmin, ymin)
  fmax <- max(xmax, ymax)
  
    fit <- lm(mergedDE$log2FoldChange.y~mergedDE$log2FoldChange.x)
  print(summary(fit))
  
  p <- ggplot(data = mergedDE, aes(x = salmon.DESeq2.FC.x, y = salmon.DESeq2.FC.y)) + 
    geom_point(size = .5, alpha = .5, col = "black") + geom_smooth(method = "lm", col = "gray") + 
    geom_vline(aes(xintercept = 0), lty = "dotted") + geom_hline(aes(yintercept = 0), lty = "dotted") +
    #geom_abline(aes(intercept = 0, slope = 1), lty = "dotted", col = "darkgrey") + 
    xlab("Acetone, 24 hours") + ylab("Phenanthrene, 24 hours") + 
    theme_bw(base_size = 12) + xlim(fmin, fmax) + ylim(fmin, fmax) + 
    #labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                   "Intercept =",signif(fit$coef[[1]],5 ),
    #                   " Slope =",signif(fit$coef[[2]], 5),
    #                   " P =",signif(summary(fit)$coef[2,4], 5))) + 
    coord_fixed()
  
  #print(summary(lm(mergedDE$log2FoldChange.x ~ mergedDE$log2FoldChange.y)))
  print(cor.test(mergedDE$salmon.DESeq2.FC.x, mergedDE$salmon.DESeq2.FC.y, method = "pearson"))
  
ggsave("Eve_corr.svg", pv, width = 8, height = 8)
ggsave("Ecy_corr.svg", pc, width = 4, height = 4)
ggsave("Gla_corr.svg", pl, width = 4, height = 4)
