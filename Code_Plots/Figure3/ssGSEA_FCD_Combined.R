library(GSVA)
library(GSEABase)
library(readxl)
library(rtracklayer)

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


# ********* ssGSEA Hallmarks genesets
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/h.all.v7.5.1.symbols.gmt")
names(h.all)

# mTORC1 pathway - all 200 genes found in the RNA-seq dataset
sum(geneIds(h.all[[26]]) %in% rownames(expr))


gs_gsva <- gsva(as.matrix(expr), h.all, method = "gsva")
gs_gsva <- t(gs_gsva)

gs_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
gs_ssgsea <- t(gs_ssgsea)

save(file="genesets_hallmarks_enrichment.Rd", 
     list = c("gs_gsva", "gs_ssgsea"))


# ******** ssGSEA Canonical pathways
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c2.cp.v7.5.1.symbols.gmt")
names(h.all)

c2_gsva <- gsva(as.matrix(expr), h.all, method = "gsva")
c2_gsva <- t(c2_gsva)


c2_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
c2_ssgsea <- t(c2_ssgsea)

save(file="genesets_Pathways_enrichment.Rd", 
     list = c("c2_gsva", "c2_ssgsea"))

# ******** ssGSEA Cell types pathways
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c8.all.v7.5.1.symbols.gmt")
names(h.all)

c8_gsva <- gsva(as.matrix(expr), h.all, method = "gsva")
c8_gsva <- t(c8_gsva)


c8_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
c8_ssgsea <- t(c8_ssgsea)

save(file="genesets_CellTypes_enrichment.Rd", 
     list = c("c8_gsva", "c8_ssgsea"))

# ******** ssGSEA Gene Ontology
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c5.go.v7.5.1.symbols.gmt")
names(h.all)

c5_gsva <- gsva(as.matrix(expr), h.all, method = "gsva")
c5_gsva <- t(c5_gsva)


c5_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
c5_ssgsea <- t(c5_ssgsea)

save(file="genesets_GOterms_enrichment.Rd", 
     list = c("c5_gsva", "c5_ssgsea"))


# ******** ssGSEA TF targets
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c3.tft.gtrd.v7.5.1.symbols.gmt")
names(h.all)

c3_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
c3_ssgsea <- t(c3_ssgsea)

save(file="genesets_TFT_enrichment.Rd", 
     list = c("c3_ssgsea"))

# ******** ssGSEA Human Phenotype Ontology (HPO)
h.all <- getGmt("/media/hdd/FCM-UNICAMP/assets/MSigDB/c5.hpo.v7.5.1.symbols.gmt")
names(h.all)

c5_ssgsea <- gsva(as.matrix(expr), h.all, method = "ssgsea")
c5_ssgsea <- t(c5_ssgsea)

save(file="genesets_HPO_enrichment.Rd", 
     list = c("c5_ssgsea"))
