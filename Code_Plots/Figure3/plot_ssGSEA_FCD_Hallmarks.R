
library(dplyr)
library(pheatmap)
library(RColorBrewer)

rm(list = ls())

load("genesets_hallmarks_enrichment.Rd")


load("FCD_Expression_BatchCorrection.rds")

head(pheno_comb)
head(gs_ssgsea)

pval <- rep(NA, ncol(gs_ssgsea))
for (i in 1:ncol(gs_ssgsea)) {
  
  s1 <- pheno_comb %>% filter(Subtype == "FCD IIb") %>% pull(sample)
  s2 <- pheno_comb %>% filter(Subtype == "FCD IIa") %>% pull(sample)
  res <- wilcox.test(gs_ssgsea[s1, i], gs_ssgsea[s2, i])
  pval[i] <- res$p.value
  
}

df_test <- data.frame(Sig = colnames(gs_ssgsea),
                      pval = pval)
df_test$pval_adj <- p.adjust(df_test$pval, method = "BH")

sig <- df_test %>% filter(pval_adj < 0.05) %>% pull(Sig)


colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

annotation_rows = data.frame(
  Subtype = factor(pheno_comb$Subtype),
  Cohort = factor(pheno_comb$Cohort)
)
rownames(annotation_rows) <- row.names(gs_ssgsea)

ann_colors = list(
  #Donor = c("SLE392" = "tan", "SLE413" = "tan1", "SLE298" = "tan3", "SLE456" = "tan4"),
  Subtype = c("FCD IIa" = "magenta", "FCD IIb" = "navyblue"),
  Cohort = c(Aus = "tan", Br = "tan3")
)

pheatmap(t(gs_ssgsea[, sig]),
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
