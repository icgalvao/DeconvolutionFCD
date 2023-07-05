rm(list = ls())


homedir <- Sys.getenv("google_drive_i199200")

setwd(paste0(homedir,'Meu Drive/Doutorado/Deconvulucao/results/Histoquimica/FCD/GFAP'))

# Load libraries 
library(tidyverse)
library(readxl)

# Read data on immunochemistry comparing FCD IIa and IIb 
table <- read_excel('GFAP.xlsx', col_types = c("text", "numeric", "text",
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", "skip"))


# Data manipulating 

str(table)

data <- table %>% 
          rename('group' = 'GRUPO', 
                 'n_measures' = 'N# de medidas',
                 'pct_area' = '% Area',
                 'min_thr' = 'MinThr',
                 'max_thr' = 'MaxThr') %>% 
          filter(!is.na(group)) %>% 
          mutate(subtype = ifelse(group == "2A", "FCD IIa", "FCD IIb"))

colnames(data) <- tolower(colnames(data))


# Plotting 

# Calculate n, mean, and standard deviation of the sample, and 
# then calculate the margin of error 
data_barplot <- data %>% 
  group_by(subtype) %>% 
  summarize(n = n(), 
            mean = mean(pct_area),
            s = sd(pct_area)) %>% 
  mutate(margin = qt(0.975, df=n-1) * s/sqrt(n),
         li = mean - margin,
         ui = mean + margin)

# plot with mean and confidence interval OF THE MEAN 
barplot <- ggplot(data_barplot, aes(x = subtype, y = mean, fill = subtype, color = subtype)) + 
  geom_col(alpha = 0.9, width = 0.7, position = "dodge") +
  geom_errorbar(aes(ymin = li, ymax = ui), width = 0.3, color = "black" ) +
  #geom_point(data = data, aes(y = pct_area), position = position_jitter(width = 0.2)) +
  labs(x = "", y = "% Area", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c("purple", "orange")) + 
  scale_color_manual(values = c("purple", "orange")) + 
  scale_y_continuous(limits = c(0, 43))  + # to draw asterisks later 
  theme_bw() + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 14, face = "plain", color = "black"),
        legend.title = element_text(size = 14, face = "plain", color = "black"),
        legend.text = element_text(size = 12, face = "plain", color = "black"), 
        panel.grid.major = element_blank()
  )


barplot

# Saving 
#ggsave(filename = "barplot_GFAP.svg", plot = barplot, dpi = 300)

# Statistical test 
# Simple wilcox test 
test <- wilcox.test(data$pct_area ~ data$group)
test

save.image("barplot_GFAP_clean.RData")
