
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)

rm(list = ls())

load("genesets_Pathways_enrichment.Rd")


load("FCD_Expression_BatchCorrection.rds")

head(pheno_comb)
head(c2_ssgsea)

pval <- rep(NA, ncol(c2_ssgsea))
for (i in 1:ncol(c2_ssgsea)) {
  
  s1 <- pheno_comb %>% filter(Subtype == "FCD IIb") %>% pull(sample)
  s2 <- pheno_comb %>% filter(Subtype == "FCD IIa") %>% pull(sample)
  res <- wilcox.test(c2_ssgsea[s1, i], c2_ssgsea[s2, i])
  pval[i] <- res$p.value
  
}

df_test <- data.frame(Sig = colnames(c2_ssgsea),
                      pval = pval)

df_test$pval_adj <- p.adjust(df_test$pval, method = "BH")
sig <- df_test %>% filter(pval_adj < 0.09,
                          str_detect(Sig, "^BIOCARTA_")) %>% pull(Sig)

sig <- df_test %>% filter(pval_adj < 0.1,
                          str_detect(Sig, "^KEGG_")) %>% pull(Sig)

sig <- df_test %>% filter(pval_adj < 0.1,
                          str_detect(Sig, "^PID_")) %>% pull(Sig)


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
rownames(annotation_rows) <- row.names(c2_ssgsea)

ann_colors = list(
  #Donor = c("SLE392" = "tan", "SLE413" = "tan1", "SLE298" = "tan3", "SLE456" = "tan4"),
  Subtype = c("FCD IIa" = "magenta", "FCD IIb" = "navyblue"),
  Cohort = c(Aus = "tan", Br = "tan3")
)

pheatmap(t(c2_ssgsea[, sig]),
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
