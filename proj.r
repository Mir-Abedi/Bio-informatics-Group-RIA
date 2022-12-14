# setRepositories(graphics = getOption("menu.graphics"), ind = NULL, addURLs = character()) run this line and select 2
# install.packages('GEOquery')
# install.packages('limma')
# install.packages("ggpubr")
# install.packages("umap")
# install.packages('Rtsne')
# install.packages("pheatmap")


library(GEOquery)
library(Biobase)
library(magrittr)
library(dplyr)
library(ggpubr)
library(umap)
library(ggplot2)
library(Rtsne)
library(limma)
library(pheatmap)



setwd("G:\\UNI\\Resources\\Fall 2022\\Bio\\Project\\Bio-informatics-Group-RIA")
data <- getGEO('GSE48558', GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = 'resources\\')

mat <- exprs(data[[1]])
dim(mat)
View(mat)


boxplot(mat)

###

groups = data[[1]]$source_name_ch1 == "AML Patient"

umap_reduction <- umap(t(mat))

ggplot(data.frame(umap_reduction$layout), aes(X1, X2, color=groups)) + geom_point()

### 

pca_resuction = prcomp(mat, scale = TRUE, center = TRUE, retx = T)

ggplot(data.frame(pca_resuction$rotation), aes(PC1, PC2, color=data[[1]]$source_name_ch1)) + geom_point()

dim(pca_resuction$rotation)
pca_resuction$rotation
data.frame(pca_resuction$rotation)['PC1']

###

tsne_reduction = Rtsne(normalize_input(t(mat)))
ggplot(data.frame(tsne_reduction$Y), aes(X1, X2, color=data[[1]]$source_name_ch1)) + geom_point()


##

gset = data[[1]]

gs = factor(groups)
different_groups = make.names(c("N", "T"))
levels(gs) = different_groups
gset$group = gs
design = model.matrix(~~group + 0, gset)
colnames(design) = levels(gs)

fit = lmFit(gset, design)
cont.matrix = makeContrasts(contrasts = "T-N", levels = design)
fit = contrasts.fit(fit, cont.matrix)
fit = eBayes(fit, 0.01)
tab = topTable(fit, adjust="fdr", sort.by = "B", number = Inf)

###

cors = cor(mat)
png("cors.png")
pheatmap(cors, labels_row = groups + 0, labels_col = groups + 0, border_color = NA, fontsize = 6)
dev.off()

