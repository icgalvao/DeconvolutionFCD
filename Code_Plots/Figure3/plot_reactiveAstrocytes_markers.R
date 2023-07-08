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


load("FCD_Expression_BatchCorrection.rds")
dim(expr_mat)
head(expr_mat)


colnames(expr_mat)
rownames(expr_mat)


# Gencode GTF with gene annotation
gtf <- import.gff("/media/hdd/FCM-UNICAMP/assets/Gencode/gencode.v40.annotation.gtf.gz")
head(gtf)
genes <- gtf[gtf$type %in% "gene", ]
genes_sort <- genes[match(rownames(expr_mat), genes$gene_id), ]

tail(genes_sort$gene_id)
tail(rownames(expr_mat))


res <- lapply(seq_along(rownames(expr_mat)), function(i){
  aux <- strsplit(rownames(expr_mat)[i], "\\.")
  aux[[1]][1]
  
})

df_exp <- as.data.frame(expr_mat)
head(df_exp)
df_exp$gene_id_simp <- unlist(res)
df_exp$gene_name <- toupper(genes_sort$gene_name)

markers <- readxl::read_xlsx(path = "Reactive_Astrocytes_Markers.xlsx") %>%
  as.data.frame()
head(markers)


head(df_exp)
head(expr_mat)

res <- df_exp %>% filter(gene_id_simp %in% markers$EnsemblID)
#res <- df_exp %>% filter(gene_name %in% c("GFAP", "VIM"))
rownames(res) <- res$gene_name
res_mat <- as.matrix(res[, 1:19])

colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

annoation_cols = data.frame(
  Subtype = factor(pheno_comb$Subtype),
  Cohort = factor(pheno_comb$Cohort)
  
)
rownames(annoation_cols) <- colnames(res_mat)


annoation_rows = data.frame(
  Astrocytes = factor(markers$Regulation),
  Function = factor(markers$Function)
)
rownames(annoation_rows) <- markers$Gene_symbol

show_col(pal_futurama("planetexpress")(12))
function_colors <- pal_futurama("planetexpress")(12)[7:12]
#function_colors <- brewer.pal(6, "Purples")
function_colors <- stepped2(n=6)
names(function_colors) <- levels(factor(markers$Function))



ann_colors = list(
  Subtype = c("FCD IIa" = "purple", "FCD IIb" = "orange"),
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

