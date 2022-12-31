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



setwd("/Users/siasor88/Documents/GitHub/Bio-informatics-Group-RIA")
data <- getGEO('GSE48558', GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = 'resources/')

mat <- exprs(data[[1]])
dim(mat)
View(mat)


boxplot(mat)

###

groups = data[[1]]$source_name_ch1 == "AML Patient"

umap_reduction <- umap(t(mat))
View(umap_reduction$layout)
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
tablee = subset(topTable(fit, adjust="fdr", sort.by = "B", number = Inf), 
             select = c("ID",
                        "adj.P.Val",
                        "P.Value",
                        "t",
                        "B",
                        "logFC",
                        "Gene.symbol",
                        "Gene.title"
                        )
             )

###

cors = cor(mat)
png("cors.png")
pheatmap(cors, labels_row = groups + 0, labels_col = groups + 0, border_color = NA, fontsize = 6)
dev.off()


#Phase 2
meaningful_treshold = 0.05


# Genes with High regulation

up_genes<- subset(tablee, (logFC > 1) & (adj.P.Val < meaningful_treshold))

up_genes$Gene.symbol   
# as you can see some of the symbols are concated by /// patter
#Seperating ///:
gene_symbols<- unique(
                    as.character(
                      strsplit2(
                        unique(up_genes$Gene.symbol),   # Symbols from data with out seperating /// data
                        "///")
                      )
                    )
# removing empty lable
gene_symbols <- gene_symbols[gene_symbols != '']

write.table(gene_symbols, file="Phase2/AML_Up_Genes.tsv", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Genes with Low regulation
down_gene<- subset(tablee, (logFC < -1) & (adj.P.Val < meaningful_treshold))

#Repeating Same process as above
down_gene$Gene.symbol   
# as you can see some of the symbols are concated by /// patter
#Seperating ///:
gene_symbols<- unique(
  as.character(
    strsplit2(
      unique(down_gene$Gene.symbol),   # Symbols from data with out seperating /// data
      "///")
  )
)
# removing empty lable
gene_symbols <- gene_symbols[gene_symbols != '']

write.table(gene_symbols, file="Phase2/AML_Down_Genes.tsv", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



