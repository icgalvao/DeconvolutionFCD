library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)

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



markers_plot <- c("GFAP", "Nestin", "Vimentin")
gids <- markers %>% filter(Marker %in% markers_plot) %>% 
  pull(EnsemblID) %>% as.vector()

df_plot <- df_exp_annot %>% filter(EnsemblID %in% gids)

df_long <- reshape2::melt(df_plot, id = c("EnsemblID", "Marker", "Function", "Gene_symbol",
                                "Regulation"), 
                variable.name = "sample")

df_long_annot <- inner_join(df_long, pheno, by = c("sample" = "ID"), keep = T)
head(df_long_annot)


ggplot(df_long_annot,
       aes(x=Gene_symbol, y=value,
           fill = PA,
           color = PA)) + 
  geom_violin(position=position_dodge(0.8), 
              trim = FALSE, alpha = 0.7, na.rm = TRUE) +
  geom_jitter(position=position_dodge(0.8)) +
  #geom_boxplot(outlier.shape = NA, na.rm = TRUE, color = "black", fill = "white", width = 0.25) +
  stat_summary(position=position_dodge(0.8), fun=mean, 
               geom="crossbar", size=0.5, width = 0.2, color="black", na.rm = TRUE) +
  #scale_y_log10(labels=trans_format('log10', math_format(10^.x))) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(face="plain", colour="black", size=12),
        axis.title.y = element_text(face="plain", colour="black", size=12),
        #legend.title=element_blank(),
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major=element_blank(),
        legend.text = element_text(face="plain", colour="black", size=12),
        legend.title = element_text(face="plain", colour="black", size=12),
        legend.position = "left",
        plot.title = element_text(face="plain", colour="black", size=14, hjust = 0)
  ) +
  scale_fill_manual(values = c("purple", "orange") ) +
  scale_color_manual(values = c("purple", "orange") ) +
  expand_limits(y = c(5, 20)) +
  #facet_wrap(~ gene_name, ncol = 4, scales = "free_y") +
  xlab("") +
  ylab("Normalized gene expression") +
  ggtitle("RNA-seq validation")


res <- lapply(gids, function(id){
  res <- wilcox.test(value ~ PA, data = df_long_annot,
                     subset = EnsemblID %in% id,
                     exact = T)
  symbol <- markers$Marker[match(id, markers$EnsemblID)]
  
  cat(id, symbol, res$p.value, "\n")
  
})



