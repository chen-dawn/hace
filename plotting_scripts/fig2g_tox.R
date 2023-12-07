library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)
library(data.table)

##### Editing is long range and directional. #####

out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"


 df <- read_csv("~/Dropbox (Harvard University)/03Helicase/CellTiterGlo/230721_Helicase_ToxD3_Redo.csv")

max_val <- 78129.91167
df <- df %>% 
  mutate(Norm1 = Norm1/max_val) %>%
  mutate(Norm2 = Norm2/max_val) %>%
  mutate(Norm3 = Norm3/max_val) 
  
df$avg <- apply(df[,c("Norm1", "Norm2", "Norm3")], 1, mean)
df$sd <- apply(df[,c("Norm1", "Norm2", "Norm3")], 1, sd)
df <- df %>% filter(!(ConditionDisplay %in% c("pUC19", "Ns3h", "PcrA", "PcrA M6", "BLM"))) %>% 
  filter(!is.na(ConditionDisplay))
  

# ggplot(df, aes(ConditionDisplay, avg, group = type, fill = type)) + geom_bar(stat = "identity")
p1 <- ggplot(df, aes(ConditionDisplay, avg, fill = type)) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=avg-sd/2, ymax=avg+sd), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = "dodge") + 
  # geom_point(data=plotting_subset_points, aes(helicaseDisplay, editRate, group=loci), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  facet_wrap(~guide, nrow = 2) + 
  geom_hline(yintercept = 1) + 
  xlab("Construct") + 
  ylab("Cell Viability") + 
  theme(axis.line = element_line(colour = "black"),
        
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = c(0.1, 0.82)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig2g_toxicity.pdf"), p1, units = "cm", width = 16, height = 14, dpi = 300)
