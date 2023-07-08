library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)

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
head(expr_mat)
range(expr_mat) #log2 normalized

head(pheno_comb)
samples_2a <- pheno_comb %>% filter(Subtype == "FCD IIa") %>%
  pull(sample) %>% as.vector()
samples_2b <- pheno_comb %>% filter(Subtype == "FCD IIb") %>%
  pull(sample) %>% as.vector()

# Create ranking - LFC 2a vs 2b

FCD_rank <- rep(NA, nrow(expr_mat))
names(FCD_rank) <- rownames(expr_mat)
for (i in 1:nrow(expr_mat)) {
  expr_a <- mean(expr_mat[i, colnames(expr_mat) %in% samples_2a], na.rm = T)
  expr_b <- mean(expr_mat[i, colnames(expr_mat) %in% samples_2b], na.rm = T)
  FCD_rank[i] <- expr_b - expr_a
}

head(FCD_rank)
FCD_rank = sort(FCD_rank, decreasing = TRUE)

markers <- readxl::read_xlsx(path = "Reactive_Astrocytes_Markers.xlsx") %>%
  as.data.frame()
head(markers)

glist <- list("Reactive_Astrocytes_Markers" = markers$EnsemblID)

glist$Reactive_Astrocytes_Markers

fgseaRes <- fgsea(pathways = glist, 
                  stats    = FCD_rank,
                  minSize  = 15,
                  maxSize  = 500)

plotEnrichment(glist[["Reactive_Astrocytes_Markers"]],
               FCD_rank) + labs(title="Reactive astrocytes markers - FCD IIb vs FCD IIa ")
