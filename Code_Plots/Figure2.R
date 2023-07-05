library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(readxl)

rm(list = ls())
# Deconvolution using Multiple sclerosis dataset
# # MS lesions and controls - scRNA-seq 10X
# # 48,920 cells
# FCD type 2a and 2b

homedir <- Sys.getenv("google_drive_i199200")
setwd(paste0(homedir, "Meu Drive/Doutorado/Deconvulucao/results/CIBERSORTx/CSX_fractions/FCD/MS_signature"))

# CibersortX Analysis with Relative results 
dat1 <- read.table(file = "AusFabioFiltered/output/CIBERSORTx_Adjusted.txt",
                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")
dat1$Dataset <- "Dataset 1"

dat2 <- read.csv(file = "Holanda/output/CIBERSORTx_Adjusted.txt", 
                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")
dat2$Dataset <- "Dataset 2"

dat <- rbind(dat1, dat2) %>% select(!c("P.value", "Correlation", "RMSE"))
head(dat)

dat <- melt(data = dat)
colnames(dat) <- c("sample", "Dataset", "pop", "fraction")
head(dat, 20)

# Pheno RNA-seq data
pheno1 <- read_excel(paste0(homedir, "Meu Drive/Doutorado/Deconvulucao/data/metadata_padronizada/std_meta_FCD_AusFabioFiltered.xlsx"))
pheno1 <- select(pheno1, c("sample", "diagnosis"))
head(pheno1)
colnames(pheno1) <- c("sample", "Subtype")

pheno2 <- read_excel(paste0(homedir, "Meu Drive/Doutorado/Deconvulucao/data/metadata_padronizada/std_meta_Holanda.xlsx")) %>% 
  filter(diagnosis %in% c("FCD IIa", "FCD IIb")) %>% 
  select(sample, diagnosis)
head(pheno2)
colnames(pheno2) <- c("sample", "Subtype")

pheno_comb <- rbind(pheno1,  pheno2)

dat_annot <- dplyr::inner_join(dat, pheno_comb, by = "sample")
head(dat_annot, 20)

dat_annot <- dat_annot %>%
  mutate(fraction_perc = fraction*100)


dat_annot <- dat_annot %>% 
  mutate(Pop_Label = case_when(
    pop == 'EN.L5.6' ~ 'EN L5/6',
    pop == 'EN.L2.3.B' ~ 'EN L2/3 B',
    pop == 'IN.PVALB' ~ 'IN PVALB',
    pop == 'EN.L4' ~ 'EN L4',
    pop == 'EN.MIX' ~ 'EN MIX',
    pop == 'IN.SST' ~ 'IN SST',
    pop == 'IN.VIP' ~ 'IN VIP',
    pop == 'EN.PYR' ~ 'EN PYR',
    pop == 'Endo.cells' ~ 'Endo.cells',
    pop == 'OPC' ~ 'OPC',
    pop == 'EN.L2.3.A' ~ 'EN L2/3 A',
    pop == 'IN.SV2C' ~ 'IN SV2C',
    pop == 'OL.A' ~ 'OL-A',
    pop == 'OL.B' ~ 'OL-B',
    pop == 'OL.C' ~ 'OL-C',
    pop == 'Astrocytes' ~ 'Astrocytes',
    pop == 'Stromal.cells' ~ 'Stromal cells',
    pop == 'Glia.MIX' ~ 'Glia-MIX',
    pop == 'Microglia' ~ 'Microglia',
    pop == 'Phagocytes' ~ 'Phagocytes',
    pop == 'T.cells' ~ 'T-cells',
    pop == 'B.cells' ~ 'B-cells'
    
    
    
  ) )

head(dat_annot)


# *** Excitatory and Inhibitory neurons
sub_dat <- dat_annot %>% filter(str_detect(Pop_Label, "^[E|I]N "))

unique(sub_dat$Pop_Label)

sub_dat$Pop_Label <- factor(sub_dat$Pop_Label, 
                        levels = c("EN L2/3 A", "EN L2/3 B", "EN L4", "EN L5/6",
                                   "EN PYR", "EN MIX", "IN PVALB", "IN SST",
                                   "IN VIP", "IN SV2C"))

ggplot(sub_dat %>% filter(Subtype %in% c("FCD IIa", "FCD IIb")),
       aes(x=Dataset, 
           y=fraction,
           fill = Subtype,
           color = Subtype)) + 
  geom_violin(position=position_dodge(0.8),
              trim = FALSE, alpha = 0.7, na.rm = TRUE) +
  geom_boxplot(position=position_dodge(0.8), 
               outlier.shape = NA, 
               na.rm = TRUE, 
               width = 0.25,
               alpha = 0.7,
               color = "black"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size=11),
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
  xlab("") +
  ylab("Proportion") +
  facet_wrap(~ Pop_Label, ncol = 4, scales = "free_y")



# Statistical by study
res <- lapply(levels(dat_annot$pop), function(apop){
  
  
  res <- wilcox.test(fraction ~ Subtype, 
                     data = dat_annot %>% filter(Dataset %in% "Dataset 1"),
                     subset = pop %in% apop, 
                     exact = T)
  
  cat(apop, "Dataset 1",  res$p.value, "\n", 
      sep = "  ")
  
  res <- wilcox.test(fraction ~ Subtype, 
                     data = dat_annot %>% filter(Dataset %in% "Dataset 2"),
                     subset = pop %in% apop, 
                     exact = T)
  
  cat(apop, "Dataset 2",  res$p.value, "\n", 
      sep = "  ")
  
})


# Non-neuronal

sub_dat <- dat_annot %>% 
  filter(!str_detect(Pop_Label, "^EN ")) %>%
  filter(!str_detect(Pop_Label, "^IN "))

unique(sub_dat$Pop_Label)

ggplot(sub_dat %>% filter(Subtype %in% c("FCD IIa", "FCD IIb")),
       aes(x=Dataset, 
           y=fraction,
           fill = Subtype,
           color = Subtype)) + 
  geom_violin(position=position_dodge(0.8),
              trim = FALSE, alpha = 0.7, na.rm = TRUE) +
  geom_boxplot(position=position_dodge(0.8), 
               outlier.shape = NA, 
               na.rm = TRUE, 
               width = 0.25,
               alpha = 0.7,
               color = "black"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size=11),
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
  xlab("") +
  ylab("Proportion") +
  facet_wrap(~ Pop_Label, ncol = 4, scales = "free_y")



