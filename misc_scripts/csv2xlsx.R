setwd("~/Documents/Paper1_stresses/data_tables/reduced_each_by_own/DE")

#install.packages("xlsx")
library(xlsx)

nfiles <- sum(grepl(".csv", dir()))

index <- 1

for (file in dir()) {
  if (grepl(".csv", file)) {
    temp <- read.csv(file)
    ifappend = index > 1
    write.xlsx(temp, "DE.xlsx", append = ifappend, sheetName = file)
    index <- index + 1
  }
}
  