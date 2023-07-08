library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)

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

rownames(expr_mat) <- unlist(res)
df_exp <- as.data.frame(expr_mat)
df_exp$gene_id <- unlist(res)
head(df_exp)
df_exp$gene_name <- toupper(genes_sort$gene_name)

markers <- readxl::read_xlsx(path = "Reactive_Astrocytes_Markers.xlsx") %>%
  as.data.frame()
head(markers)

markers_plot <- c("GFAP", "Nestin", "Vimentin")
gids <- markers %>% filter(Marker %in% markers_plot) %>% 
  pull(EnsemblID) %>% as.vector()

df_plot <- df_exp %>% filter(gene_id %in% gids)

df_long <- melt(df_plot, id = c("gene_name", "gene_id"), 
                variable.name = "sample")

df_long_annot <- inner_join(df_long, pheno_comb, by = "sample", keep = T)
head(df_long_annot)


ggplot(df_long_annot,
       aes(x=gene_name, y=value,
           fill = Subtype,
           color = Subtype)) + 
  geom_violin(position=position_dodge(0.6), 
              trim = FALSE, alpha = 0.7, na.rm = TRUE) +
  geom_jitter(position=position_dodge(0.6)) +
  #geom_boxplot(outlier.shape = NA, na.rm = TRUE, color = "black", fill = "white", width = 0.25) +
  stat_summary(position=position_dodge(0.6), fun=mean, 
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
  ggtitle("RNA-seq discovery")



res <- lapply(gids, function(id){
  res <- wilcox.test(value ~ Subtype, data = df_long_annot,
                     subset = gene_id %in% id,
                     exact = T)
  symbol <- markers$Marker[match(id, markers$EnsemblID)]
  
  cat(id, symbol, res$p.value, "\n")
  
})
