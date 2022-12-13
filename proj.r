# setRepositories(graphics = getOption("menu.graphics"), ind = NULL, addURLs = character()) run this line and select 2
# install.packages('GEOquery')
# install.packages('limma')
# install.packages("ggpubr")


library(GEOquery)
library(Biobase)
library(magrittr)
library(dplyr)
library(ggpubr)



setwd("G:\\UNI\\Resources\\Fall 2022\\Bio\\Project\\Bio-informatics-Group-RIA")
data <- getGEO('GSE48558', GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = 'resources\\')

mat <- exprs(data[[1]])
dim(mat)
View(mat)


boxplot(mat)

###

groups = data[[1]]$source_name_ch1 == "AML Patient"

mds <- mat %>%dist()%>%cmdscale() %>% as_tibble()

colnames(mds) <- c("PC1", "PC2")
# Plot MDS
ggscatter(mds, x = "PC1", y = "PC2", 
          label = groups(),
          size = 1,
          repel = TRUE)