
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)

rm(list = ls())

load("genesets_GOterms_enrichment.Rd")


load("FCD_Expression_BatchCorrection.rds")

head(pheno_comb)
head(c5_ssgsea)

c5_ssgsea <- t(c5_ssgsea)
pval <- rep(NA, ncol(c5_ssgsea))
for (i in 1:ncol(c5_ssgsea)) {
  
  s1 <- pheno_comb %>% filter(Subtype == "FCD IIb") %>% pull(sample)
  s2 <- pheno_comb %>% filter(Subtype == "FCD IIa") %>% pull(sample)
  res <- wilcox.test(c5_ssgsea[s1, i], c5_ssgsea[s2, i])
  pval[i] <- res$p.value
  
}

df_test <- data.frame(Sig = colnames(c5_ssgsea),
                      pval = pval)

df_test$pval_adj <- p.adjust(df_test$pval, method = "BH")

sig <- df_test %>% filter(pval_adj < 0.05) %>% pull(Sig)


# BP
sig <- df_test %>% filter(pval_adj < 0.05,
                          str_detect(Sig, "^GOBP_")) %>% pull(Sig)

# CC
sig <- df_test %>% filter(pval_adj < 0.05,
                          str_detect(Sig, "^GOCC_")) %>% pull(Sig)

# MF
sig <- df_test %>% filter(pval_adj < 0.05,
                          str_detect(Sig, "^GOMF_")) %>% pull(Sig)


# sig <- df_test %>% filter(pval_adj < 0.07,
#                           str_detect(Sig, "^REACTOME_")) %>% pull(Sig)
# 
# sig <- df_test %>% filter(pval_adj < 0.07,
#                           str_detect(Sig, "^WP_")) %>% pull(Sig)



colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

annotation_rows = data.frame(
  Subtype = factor(pheno_comb$Subtype),
  Cohort = factor(pheno_comb$Cohort)
)
rownames(annotation_rows) <- row.names(c5_ssgsea)

ann_colors = list(
  #Donor = c("SLE392" = "tan", "SLE413" = "tan1", "SLE298" = "tan3", "SLE456" = "tan4"),
  Subtype = c("FCD IIa" = "magenta", "FCD IIb" = "navyblue"),
  Cohort = c(Aus = "tan", Br = "tan3")
)

pheatmap(t(c5_ssgsea[, sig]),
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         #annotation_row = annotation_rows,
         annotation_col = annotation_rows,
         show_colnames = F,
         annotation_colors = ann_colors,
         col=colors,
         scale="row",
         clustering_method = "mcquitty",
         fontsize = 9,
         fontsize_row = 9
)
