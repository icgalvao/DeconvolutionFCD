library(rtracklayer)
library(dplyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(scales)
library(pals)
library(ggsci)

rm(list = ls())

pheno <- read.table("./RNA-seq_FCD_Netherlands/annotation.txt", stringsAsFactors = FALSE,
                    header = TRUE, sep = "\t")
head(pheno)
pheno <- pheno %>% filter(PA %in% c("FCD_2A", "FCD_2B"))

#Reads were aligned directly to the human GRCh38 reference transcriptome (Gencode version 33) (2) 
# using Salmon v0.11.3 (3). Transcript counts were summarized to the gene level and scaled used 
# library size and average transcript length using the R package tximport (4). 
# Genes detected in less than 20% of the samples in any diagnosis and with counts less than 
# 6 across all samples were filtered out. The gene counts were than normalized using the weighted 
# trimmed mean of M-values (TMM) method using the R package edgeR (5). 
# The normalized counts were than log2 transformed using the voom function from the R package limma (6).

# TMM normalized + log2
expr_mat <- read.table("RNA-seq_FCD_Netherlands/expression_matrix.txt", stringsAsFactors = FALSE, header = TRUE)
colnames(expr_mat) <- gsub("^X", "", colnames(expr_mat))
colnames(expr_mat) <- gsub("\\.", "-", colnames(expr_mat))
colnames(expr_mat) %in% pheno$ID
expr_mat <- expr_mat[, colnames(expr_mat) %in% pheno$ID]
head(expr_mat)


markers <- readxl::read_xlsx(path = "Reactive_Astrocytes_Markers.xlsx") %>%
  as.data.frame()
head(markers)


res <- expr_mat[rownames(expr_mat) %in% markers$EnsemblID, ]
df_exp <- as.data.frame(res)
df_exp$EnsemblID <- rownames(res)
head(df_exp)


df_exp_annot <- df_exp %>% inner_join(markers, df_exp, by = "EnsemblID")
head(df_exp_annot)



res_mat <- as.matrix(df_exp_annot[, 1:32])
rownames(res_mat) <- df_exp_annot$Gene_symbol

colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

annoation_cols = data.frame(
  Subtype = factor(pheno$PA)
)
rownames(annoation_cols) <- colnames(res_mat)


annoation_rows = data.frame(
  Astrocytes = factor(markers$Regulation),
  Function = factor(markers$Function)
)
rownames(annoation_rows) <- markers$Gene_symbol

function_colors <- stepped2(n=6)
names(function_colors) <- levels(factor(markers$Function))



ann_colors = list(
  Subtype = c("FCD_2A" = "purple", "FCD_2B" = "orange"),
  Cohort = c(Aus = "mistyrose1", Br = "mistyrose3"),
  Astrocytes = c(Up = "#b2182b", Down = "#2166ac"),
  Function = function_colors
)

pheatmap(res_mat,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         cluster_cols = T,
         annotation_row = annoation_rows,
         annotation_col = annoation_cols,
         show_colnames = F,
         annotation_colors = ann_colors,
         color = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
         scale="row",
         treeheight_row = 25,
         treeheight_col = 25,
         clustering_method = "mcquitty"
)
