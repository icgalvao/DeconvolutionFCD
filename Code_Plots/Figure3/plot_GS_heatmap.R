library(GSVA)
library(GSEABase)
library(readxl)
library(rtracklayer)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

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

# Remove duplicates based on gene names
df_nodup <- df_exp[!duplicated(df_exp$gene_name), ]

# Expression matrix
colnames(df_nodup)
expr <- df_nodup[ , c(-20, -21)]
rownames(expr) <- df_nodup$gene_name


h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c2.cp.v7.5.1.symbols.gmt")
names(h.all)

# gs <- h.all[[grep("REACTOME_LOSS_OF_FUNCTION_OF_MECP2_IN_RETT_SYNDROME", names(h.all))]]
gs <- h.all[[grep("BIOCARTA_LYMPHOCYTE_PATHWAY", names(h.all))]]
gs <- h.all[[grep("WP_COMPLEMENT_ACTIVATION", names(h.all))]]
gs <- h.all[[grep("REACTOME_MTOR_SIGNALLING", names(h.all))]]


h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c8.all.v7.5.1.symbols.gmt")
gs <- h.all[[grep("ZHONG_PFC_C1_ASTROCYTE", names(h.all))]]
gs <- h.all[[grep("ZHONG_PFC_MAJOR_TYPES_ASTROCYTES", names(h.all))]]
gs <- h.all[[grep("ZHONG_PFC_MAJOR_TYPES_OPC", names(h.all))]]


geneIds(gs)

expr_gs <- expr[rownames(expr) %in% geneIds(gs), ]

colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

annotation_rows = data.frame(
  Subtype = factor(pheno_comb$Subtype),
  Cohort = factor(pheno_comb$Cohort)
)
rownames(annotation_rows) <- pheno_comb$sample

ann_colors = list(
  #Donor = c("SLE392" = "tan", "SLE413" = "tan1", "SLE298" = "tan3", "SLE456" = "tan4"),
  Subtype = c("FCD IIa" = "magenta", "FCD IIb" = "navyblue"),
  Cohort = c(Aus = "tan", Br = "tan3")
)

pheatmap(expr_gs,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         #annotation_row = annotation_rows,
         annotation_col = annotation_rows,
         show_colnames = F,
         annotation_colors = ann_colors,
         col=colors,
         scale="row",
         clustering_method = "mcquitty"
)
