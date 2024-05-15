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
#### Helper Functions ####
get_revcomp_base <- function(x){
  if (x == "A"){
    return("T")
  }
  if (x == "T"){
    return("A")
  }
  if (x == "G"){
    return("C")
  }
  if (x == "C"){
    return("G")
  }
}
flip_DNA <- function(MEK1_pileup){
  MEK1_pileup$ref_base <- sapply(MEK1_pileup$ref_base, get_revcomp_base)
  Aold <- MEK1_pileup$A
  Cold <- MEK1_pileup$C
  MEK1_pileup$A <- MEK1_pileup$`T`
  MEK1_pileup$C <- MEK1_pileup$G
  MEK1_pileup$G <- Cold
  MEK1_pileup$`T` <- Aold
  MEK1_pileup$n_base <- -MEK1_pileup$n_base
  return(MEK1_pileup)
}
 
MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/"
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename)

all_files_df <- all_files_df %>%
  mutate(condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition)) %>%
  mutate(condition = gsub("_", "-", condition))

fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")

# c <- reconcat_all %>% group_by(sample) %>% summarise(avg = mean(base_sum))

all_files_df <- fread("~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")
MEK1_pileup <- all_files_df %>% 
  mutate(base_sum = (A + C + `T` + G )) 
#  filter(base_sum > 10000)
# MEK1_pileup <- MEK1_pileup %>% mutate(permname = paste0(sample, chr, n_base))

# Nick Locations:
# DC203: 48027 reverse
# RUNX1: 5115 forward
# DC329: TNF: 2478 forward 
# DC330 IL6: 1534 forward
# DC331: 15562 reverse
## More Locations
# DC573 44231 reverse
# DC576 2972 forward
# DC577 7846 forward
DC203 <- MEK1_pileup %>% filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 48027) %>% filter(chr == "10_MAP2K1") %>% mutate(loci = "MAP2K1")

RUNX1 <- MEK1_pileup %>% filter(grepl("RUNX1", sample)) %>% 
  mutate(n_base = n_base - 5115) %>% filter(chr == "RUNX1_8kb_extraction") %>% mutate(loci = "RUNX1")
RUNX1 <- flip_DNA(RUNX1)

DC329 <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 2478) %>% filter(chr == "chr6_(NC_000006)_TNF_long_extraction") %>% mutate(loci = "TNF")
DC329 <- flip_DNA(DC329)

DC330 <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 1534) %>% filter(chr == "NC_000007_-_IL6_gene") %>% mutate(loci = "IL6")
DC330 <- flip_DNA(DC330)

DC331 <- MEK1_pileup %>% filter(grepl("DC331", sample)) %>% 
  mutate(n_base = n_base - 15562) %>% filter(chr =="NC_000014_-_AKT1_gene") %>% mutate(loci = "AKT1")

DC573 <- MEK1_pileup %>% filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 44231) %>% filter(chr == "DC573_02_HEK3") %>% mutate(loci = "HEK3")

DC576 <- MEK1_pileup %>% filter(grepl("DC576", sample)) %>% 
  mutate(n_base = n_base - 2972) %>% filter(chr == "DC576_05_VEGFA") %>% mutate(loci = "VEGFA")
DC576 <- flip_DNA(DC576)

DC577 <- MEK1_pileup %>% filter(grepl("DC577", sample)) %>% 
  mutate(n_base = n_base - 7846) %>% filter(chr == "DC577_06_CD209") %>% mutate(loci = "CD209")
DC577 <- flip_DNA(DC577)

# Need to change names.
reconcat_all <- rbind(RUNX1, DC203, DC329, DC330, DC331, DC573, DC576, DC577) %>% filter(n_base > -2000 & n_base < 2000)
# with_names_edited <- reconcat_all %>% mutate(temp_name = str_extract(sample, "^.*(?=_\\d_S\\d+)")) %>%
#   mutate(temp_extension = ifelse(is.na(temp_name),NA, str_extract(sample, "_\\d_S\\d+"))) %>% 
#   mutate(temp_name = gsub("_", "-",temp_name)) %>% 
#   mutate(condition = ifelse(is.na(temp_name), condition, temp_name)) %>% 
#   select(-temp_name, -temp_extension)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all %>% 
  mutate(helicase = str_extract(condition, "DC333|DC334|DC335|pBO101|GW109|GW111|GW113")) %>% 
  mutate(Cas9 = str_extract(condition, "DC332|DC324|GW134|GW133|DC327")) %>% 
  mutate(helicaseDisplay = helicase) %>% 
  mutate(Cas9Display = Cas9)

with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC333", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW109", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC334", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW111", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC335", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW113", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "pBO101", "PcrA M6")

with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC332", "nCas9 D10A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "GW133", "nCas9 D10A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC324", "nCas9 H840A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC327", "nCas9 H840A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "GW134", "dCas9")

with_mutation_rate <- with_plasmid_annotation %>%  
  filter((A+C+G+`T`)>10000)%>%
  mutate(GtoA = A/(A+G)*100) %>% 
  mutate(CtoT = `T`/(C + `T`)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) %>% 
  # Add the control samples. 
  mutate(isControl = ifelse(condition %in% c("RUNX1-Ctrl", "DC203-Control", "DC329-Ctrl", "DC330-Ctrl", "DC331-Ctrl", "DC573-DC573", "DC576-DC573", "DC577-DC573"), 1, 0))
  
control_rates <- with_mutation_rate %>% filter(isControl  ==  1) %>% 
  group_by(chr, n_base, ref_base)%>%
  summarise(meanCtoT_ctrl = mean(CtoT), 
            meanGtoA_ctrl = mean(GtoA))

with_mut_and_control <- merge(with_mutation_rate %>% filter(isControl == 0), control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)

with_mut_background_subtracted <- with_mut_and_control %>% 
  mutate(CtoT_bg_subtracted = CtoT - meanCtoT_ctrl) %>% 
  mutate(GtoA_bg_subtracted = GtoA - meanGtoA_ctrl) %>%
  mutate(CtoA = A/(C + A)) %>%
  mutate(TtoA = A/(A + `T`))

ggplot(with_mut_background_subtracted %>% filter(ref_base.x == "G") %>% filter(Cas9Display == "dCas9"), aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) +
  ylim(0,1) +
 geom_point() + facet_grid(helicaseDisplay+Cas9Display ~ loci) + xlim(-1000,1000)

avg_mut_per_base <- with_mut_background_subtracted %>% 
  pivot_longer(cols = contains("bg_subtracted"), names_to = "editMode") %>% 
  group_by(chr, n_base, ref_base.x, editMode, condition, loci, helicaseDisplay, Cas9Display) %>% 
  summarise(editRate = mean(value), sd = sd(value)) %>%
  mutate(editMode = factor(editMode, levels = c("GtoA_bg_subtracted", "CtoT_bg_subtracted")))

# Remove the background subtraction
avg_mut_per_base <- with_mut_and_control %>% 
  mutate(CtoT_bg_subtracted = CtoT) %>% 
  mutate(GtoA_bg_subtracted = GtoA) %>%
  pivot_longer(cols = contains("bg_subtracted"), names_to = "editMode") %>% 
  group_by(chr, n_base, ref_base.x, editMode, condition, loci, helicaseDisplay, Cas9Display) %>% 
  summarise(editRate = mean(value), sd = sd(value)) %>%
  mutate(editMode = factor(editMode, levels = c("GtoA_bg_subtracted", "CtoT_bg_subtracted")))

avg_mut_per_base_filtered <- avg_mut_per_base %>% 
  filter(ref_base.x %in% c("C", "G")) 
  
# Filter out bases with too high edit rate (>5%). 
positive_sample <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-pBO101-DC332")) %>%
  filter(editRate > 5)

# Filter out bases with high background.
high_in_ctrl <- control_rates %>% filter(meanCtoT_ctrl > 0.04 | meanGtoA_ctrl > 0.04) %>% filter(chr == "DC573_02_HEK3")

##### PLOT NOW. #####
newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")
salmon_on_ice <- c("#2185C5", "#7ECEFD", "#FFF6E5", "#FF7F66","#3E454C")
salmon_on_ice2 <- c("#FF7F66","#2185C5", "#7ECEFD", "#FFF6E5","#3E454C")

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
p1a_range <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 1) + 
  ylim(0,2) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))
p1a_range

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-pBO101")) %>%
  filter(!(n_base %in% positive_sample$n_base))
# plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW113-GW133"))
p1a_range2 <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 1) + 
  ylim(0,2) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) +
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))
  
  # facet_zoom(ylim = c(0, 0.5), xlim = c(-1000,1000), zoom.size = 0.5, horizontal = FALSE, show.area = F)
g <- arrangeGrob(p1a_range, p1a_range2, ncol = 1)

# ggsave(file.path(out_dir, "fig1a_range2.pdf"), g, units = "cm", width = 12, height = 14, dpi = 300)



### Also plot for the other helicases. (THIS IS HEK3)
plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC333-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
BLM <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC334-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
Ns3h <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC335-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA_M6 <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))


g <- arrangeGrob(BLM, Ns3h, PcrA, PcrA_M6, ncol = 1)
ggsave(file.path(out_dir, "fig2g_D10A_long.pdf"), g, units = "cm", width = 7, height = 14, dpi = 300)


plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC333-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
BLM <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC334-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
Ns3h <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC334-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-pBO101-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA_M6 <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))


g <- arrangeGrob(BLM, Ns3h, PcrA, PcrA_M6, ncol = 1)
ggsave(file.path(out_dir, "fig2g_H840A_long.pdf"), g, units = "cm", width = 7, height = 14, dpi = 300)


# If we were to try just merge all and faceting:
plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC333-DC332", "DC573-DC334-DC332", "DC573-DC335-DC332", "DC573-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC332_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_D10A_long_facet.pdf"), DC332_facet, units = "cm", width = 7, height = 14, dpi = 300)

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC573-DC333-DC324", "DC573-DC334-DC324", "DC573-DC335-DC324", "DC573-pBO101-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC324_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_H840A_long_facet.pdf"), DC324_facet, units = "cm", width = 7, height = 14, dpi = 300)

g <- arrangeGrob(DC332_facet, DC324_facet, ncol = 2)
ggsave(file.path(out_dir, "fig2g_both_no_header.pdf"), g, units = "cm", width = 14.5, height = 14, dpi = 300)



### Also plot for the other helicases. (THIS IS ANOTHER LOCI)
plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW109-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
BLM <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW111-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
Ns3h <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW113-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA_M6 <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,5) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))


g <- arrangeGrob(BLM, Ns3h, PcrA, PcrA_M6, ncol = 1)
ggsave(file.path(out_dir, "fig2g_D10A_long.pdf"), g, units = "cm", width = 7, height = 14, dpi = 300)


plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW109-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
BLM <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW111-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
Ns3h <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW113-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-pBO101-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
PcrA_M6 <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))


g <- arrangeGrob(BLM, Ns3h, PcrA, PcrA_M6, ncol = 1)
ggsave(file.path(out_dir, "fig2g_H840A_long.pdf"), g, units = "cm", width = 7, height = 14, dpi = 300)


# If we were to try just merge all and faceting:
plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW109-GW133", "DC329-GW111-GW133", "DC329-GW113-GW133", "DC329-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC332_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_D10A_long_facet_TNF.pdf"), DC332_facet, units = "cm", width = 7, height = 14, dpi = 300)

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW109-DC324", "DC329-GW111-DC324", "DC329-GW113-DC324", "DC329-pBO101-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC324_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_H840A_long_facet_TNF.pdf"), DC324_facet, units = "cm", width = 7, height = 14, dpi = 300)

g <- arrangeGrob(DC332_facet, DC324_facet, ncol = 2)
ggsave(file.path(out_dir, "fig2g_both_no_header_TNF.pdf"), g, units = "cm", width = 14.5, height = 14, dpi = 300)


# 

# If we were to try just merge all and faceting:
plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC330-GW109-GW133", "DC330-GW111-GW133", "DC330-GW113-GW133", "DC330-pBO101-DC332")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC332_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_D10A_long_facet_IL6.pdf"), DC332_facet, units = "cm", width = 7, height = 14, dpi = 300)

plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC330-GW109-DC324", "DC330-GW111-DC324", "DC330-GW113-DC324", "DC330-pBO101-DC324")) %>%
  filter(!(n_base %in% positive_sample$n_base)) 
DC324_facet <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 0.8) + 
  ylim(0,3) +
  scale_x_continuous(breaks = round(seq(-1500, 1500, by = 250),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_manual(values = salmon_on_ice2) + 
  theme(axis.line = element_line(colour = "black", size = unit(0.4, 'pt')),
        axis.ticks = element_line(colour = "black", size = unit(0.4, 'pt')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87)) + facet_wrap(~condition, ncol = 1) + 
  theme(strip.text.x = element_blank())

ggsave(file.path(out_dir, "fig2g_H840A_long_facet_IL6.pdf"), DC324_facet, units = "cm", width = 7, height = 14, dpi = 300)

g <- arrangeGrob(DC332_facet, DC324_facet, ncol = 2)
ggsave(file.path(out_dir, "fig2g_both_no_header_IL6.pdf"), g, units = "cm", width = 14.5, height = 14, dpi = 300)



