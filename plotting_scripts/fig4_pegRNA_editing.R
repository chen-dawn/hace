library(tidyverse)
library(gridExtra)
library(grid)
library(gt)
library(matrixTests)
library(ggrepel)
library(gridExtra)
library("viridis")
library(ggsci)
library(ggpubr)
library(cowplot)

out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"

editing_df <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt20_SF3B1_BE_PE_Validation/SF3B1_PE_all_samples_aggregate_table.csv")
df_filtered <- editing_df %>% mutate(avgEditingRate = (EditingRateR1 + EditingRateR2 + EditingRateR3)/3) %>% 
  mutate(normalizedGFPMcherryRatio = gfp_mcherry_ratio/0.308359399) %>% 
  mutate(ratioSD = sd_ratio/0.308359399) 



ggplot(df_filtered %>% filter(avgEditingRate != 0 | condition == "01JS4"), aes(avgEditingRate, normalizedGFPMcherryRatio, label = condition)) + geom_point() + geom_label()


df_filtered_for_plot <- df_filtered %>% filter(grepl("1851|1849|01JS4|1993|1868|1682|1996", condition)) %>% filter(!grepl("DC819|DC807|DC816|DC818", condition))
g3<- ggplot(df_filtered_for_plot, aes(avgEditingRate, normalizedGFPMcherryRatio, label = condition)) + 
  geom_pointrange(aes(ymin = normalizedGFPMcherryRatio -ratioSD/2, ymax = normalizedGFPMcherryRatio + ratioSD/2)) + 
  geom_point()+ #geom_label()
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Editing Rate (%)") + 
  ylab("Reporter Fold Change") 
ggsave(file.path(out_dir, "supp5b_pegRNA_editing_foldchange.pdf"), g3, units = "cm", width = 12, height = 11, dpi = 300)

write_csv(df_filtered_for_plot, file.path(out_dir, "supp5b_pegRNA_editing_foldchange.csv"))


BE_PE_both <- read_csv("~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/supp5b_BEandPE_foldChangeEditing.csv") %>% 
  mutate(avgEditingRate = (EditingRateR1 + EditingRateR2 + EditingRateR3)/3) %>% 
  mutate(EditingType = factor(EditingType, levels = c("BE", "PE", "Control")))
subset <- BE_PE_both %>% 
  filter(!grepl("DC763", condition)) %>% 
  filter(!grepl("DC803", condition)) %>% 
  filter(!grepl("DC804", condition)) %>% 
  filter(!grepl("DC808", condition)) %>% 
  filter(!grepl("DC809", condition)) %>% 
  filter(!grepl("DC810", condition)) %>%
  filter(!grepl("DC820", condition)) 
  

m1 <- lm(normalizedGFPMcherryRatio~avgEditingRate, data = subset)
g4 <- ggplot(subset, aes(avgEditingRate, normalizedGFPMcherryRatio, label = condition, color = EditingType)) + 
  geom_pointrange(aes(ymin = normalizedGFPMcherryRatio -ratioSD/2, ymax = normalizedGFPMcherryRatio + ratioSD/2)) + 
  geom_point()+ 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_npg() + xlim(0,40) + ylim(0,15) + 
  geom_abline(slope = coef(m1)[["avgEditingRate"]], 
              intercept = coef(m1)[["(Intercept)"]])+
  # scale_color_manual(values = newPalette400_rotate) + 
  xlab("Editing Rate (%)") + 
  ylab("Reporter Fold Change") 
cor(subset$avgEditingRate, subset$normalizedGFPMcherryRatio)
# [1] 0.880843

ggsave(file.path(out_dir, "supp5b_pegRNA_editing_foldchange_colored.pdf"), g4, units = "cm", width = 12, height = 8, dpi = 300)

