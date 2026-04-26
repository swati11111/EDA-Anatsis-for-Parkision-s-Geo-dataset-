#Install & Load Packages#

if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "pheatmap", "ggplot2"))

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)

# Dwonload dataset 

gse <- getGEO("GSE49036", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])

# View dimensions
dim(exprSet)
head(exprSet[,1:5])
pheno <- pData(gse[[1]])
head(pheno)

# Check group/condition column
colnames(pheno)

boxplot(exprSet, outline = FALSE, las = 2,
        main = "Boxplot Before Normalization")

# If needed normalize
exprSet_norm <- normalizeBetweenArrays(exprSet)

boxplot(exprSet_norm, outline = FALSE, las = 2,
        main = "Boxplot After Normalization")
#View PCA PLot

pca <- prcomp(t(exprSet_norm), scale. = TRUE)

pca_df <- data.frame(pca$x, group = pheno$characteristics_ch1)

ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA Plot")

# Select top variable genes
var_genes <- apply(exprSet_norm, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]

#view Heat Map
pheatmap(exprSet_norm[top_genes, ],
         scale = "row",
         show_rownames = FALSE,
         main = "Top Variable Genes Heatmap")
plotDensities(exprSet_norm, main = "Density Plot")
sum(is.na(exprSet_norm))
summary(exprSet_norm)