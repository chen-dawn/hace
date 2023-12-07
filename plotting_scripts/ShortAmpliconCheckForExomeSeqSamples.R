library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)
library(data.table)
##### Quantifying editing rates. #####

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

MEK1_dir <- "/Volumes/broad_thechenlab/NextSeqOutput_Xerneas/230703_VH01286_19_AACJG2NM5/Data/Intensities/BaseCalls/exome_24_well"

#### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename) %>%
  mutate(base_sum = (A + C + `T` + G )) %>%
  filter(base_sum > 100)

MEK1_pileup <- all_files_df
MEK1_pileup <- MEK1_pileup %>% 
  mutate(base_sum = (A + C + `T` + G )) %>% 
  filter(base_sum > 10000)

with_mutation_rate <- MEK1_pileup %>%  
  filter((A+C+G+`T`)>10000)%>%
  mutate(GtoA = A/(A+G)*100) %>% 
  mutate(CtoT = `T`/(C + `T`)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) 

ggplot(with_mutation_rate %>% filter(ref_base == "G"), aes(n_base, GtoA)) + ylim(0,2) +
  geom_point() + 
  facet_wrap(~ condition, nrow = 2)

# Nick Locations:
# DC203: 140 reverse
# RUNX1_300bp: 146 reverse
# DC329: TNF: 217 forward 
# DC330 IL6: 215 forward
## More Locations
# DC572 130 reverse
# DC573 119 reverse
# DC576 251 forward
# DC577 263 forward
DC203 <- MEK1_pileup %>% 
  mutate(n_base = n_base - 140) %>% filter(chr == "10_MAP2K1_PCR_Product") %>% mutate(loci = "MAP2K1")


# Need to change names.
# Remove MAP2K1, seems to be in the wrong sgRNA set. 
reconcat_all <- rbind(DC203) %>% 
  filter(n_base > -150 & n_base < 150) %>% 
  mutate(condition = gsub("24well-", "", condition)) %>% 
  mutate(condition = gsub("GW109", "DC333", condition)) 
# filter(!grepl("GW109", condition))

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all %>% 
  mutate(helicase = str_extract(condition, "DC333|DC334|DC335|pBO101|DC348|GW109|GW111|GW113")) %>% 
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
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC348", "PcrA M6")

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
  mutate(isControl = ifelse(is.na(helicaseDisplay), 1, 0)) %>% 
  mutate(isBaseline = ifelse(is.na(helicaseDisplay) & is.na(Cas9Display), 1, 0))

# control_rates <- with_mutation_rate %>% filter(isBaseline  ==  1) %>% 
#   group_by(chr, n_base)%>%
#   summarise(meanCtoT_ctrl = mean(CtoT), 
#             meanGtoA_ctrl = mean(GtoA))
# 
# with_mut_and_control <- merge(with_mutation_rate, control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)
# 
# with_mut_background_subtracted <- with_mut_and_control %>% 
#   mutate(CtoT_bg_subtracted = CtoT - meanCtoT_ctrl) %>% 
#   mutate(GtoA_bg_subtracted = GtoA - meanGtoA_ctrl) %>%
#   mutate(CtoA = A/(C + A)) %>%
#   mutate(TtoA = A/(A + `T`))



# For each loci, we want to quantify the editing rate.
allC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample,condition, loci, helicaseDisplay, Cas9Display) %>% 
  dplyr::summarise(meanCtoT_bigger_than_10bp = mean(CtoT))

allC_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample) %>%
  dplyr::summarise(meanCtoT_spacer = mean(CtoT))

allC_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample) %>% 
  dplyr::summarise(meanCtoT_smaller_than_10bp = mean(CtoT))

allG_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  dplyr::summarise(meanGtoA_spacer = mean(GtoA))

allG_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  dplyr::summarise(meanGtoA_bigger_than_10bp = mean(GtoA))

allG_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  dplyr::summarise(meanGtoA_smaller_than_10bp = mean(GtoA))


ggplot(with_mutation_rate, aes(n_base, GtoA)) + geom_point() + 
  facet_wrap(~condition, nrow = 2) + 
  ylim(0,1)


# Merge all the stats
merged_stats <- list(allC_10bp, allC_bigger_than_10bp, allC_smaller_than_minus10bp, allG_10bp, allG_bigger_than_10bp, allG_smaller_than_minus10bp) %>% 
  reduce(full_join, by = "sample")

merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
  mutate(baseMode = str_extract(mode, "GtoA|CtoT")) 


# Take average across replicates.
avg_across_reps <- merged_stats_pivot_long %>% 
  filter(!is.na(editRate)) %>%
  group_by(condition, loci, helicaseDisplay, Cas9Display, baseMode, mode) %>% 
  summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) 

##### PLOT NOW. #####
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("MAP2K1")) %>% 
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(mode == "bigger_than_10bp") %>% 
  filter(baseMode == "GtoA")
# filter(!grepl("GW109", condition)) %>% 
#mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))

plotting_subset_points <- merged_stats_pivot_long %>% 
  filter(loci %in% c("MAP2K1")) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  # filter(!is.na(editRate)) %>%
  filter(mode == "bigger_than_10bp") %>% 
  filter(baseMode == "GtoA")

write_csv(plotting_subset_points, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_ShortAmpliconSeqForExomeCheck.csv")

newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")

g1 <- ggplot(plotting_subset, aes(loci, meanEditRate, fill = helicaseDisplay, group = helicaseDisplay)) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), color = "black", size = 0.4) + 
  geom_point(data=plotting_subset_points, aes(loci, editRate, group=helicaseDisplay), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_brewer(direction = -1) + 
  xlab("Loci") + 
  ylab("Mean Edits Per kb") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.9, 0.8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
