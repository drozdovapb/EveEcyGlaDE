options(stringsAsFactors = F)
library(ggplot2)

aliper <- read.csv("/homes/brauerei/polina/Documents/Paper1_stresses/alignment_percent_to_ani.csv")


aliper$species <- substr(aliper$Sample, 1, 3)



ggplot(aliper, aes(x = species, y = X * 100, fill=species)) + geom_boxplot() + 
  facet_grid(.~Reference) + ylab("Mapped raw reads, %") + xlab("Source read species") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00"))+ theme_bw(base_size = 14)
ggsave("~/Documents/Paper1_stresses/alignment_percent.svg")
