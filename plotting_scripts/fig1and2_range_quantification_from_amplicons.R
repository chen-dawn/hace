library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)
library(ggbeeswarm)
library(data.table)
##### Quantifying editing rates. #####

theme_set(theme_bw())
theme_update(axis.line = element_line(colour = "black"))
theme_update(panel.grid.major = element_blank())
theme_update(panel.grid.minor = element_blank())
theme_update(panel.border = element_blank())
theme_update(panel.background = element_blank())
theme_update(text = element_text(family = "Helvetica"))
theme_update(text = element_text(size = 9))
newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")
salmon_on_ice <- c("#2185C5", "#7ECEFD", "#FFF6E5", "#FF7F66","#3E454C")


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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt01_Multiple_locations_all_helicase_Cas9/"

##### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename)

fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt01_Multiple_locations_all_helicase_Cas9/merged_all_amplicon_pileups.csv")

MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt01_Multiple_locations_all_helicase_Cas9/merged_all_amplicon_pileups.csv")
MEK1_pileup <- MEK1_pileup %>% 
  mutate(base_sum = (A + C + `T` + G )) %>% 
  filter(base_sum > 10000)

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
DC203 <- MEK1_pileup %>% filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 140) %>% filter(chr == "10_MAP2K1_PCR_Product") %>% mutate(loci = "MAP2K1")
DC203_ctrl <- MEK1_pileup %>% filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 102) %>% filter(chr == "MEK1i1_300bp") %>% mutate(loci = "MAP2K1") %>% 
  filter(condition %in% c("DC203", "DC203-DC324", "DC203-DC332", "DC203-GW134")) %>% 
  mutate(chr = "10_MAP2K1_PCR_Product")

RUNX1 <- MEK1_pileup %>% filter(grepl("RUNX1", sample)) %>% 
  mutate(n_base = n_base - 146) %>% filter(chr == "RUNX1_300bp") %>% mutate(loci = "RUNX1")

DC329 <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 217) %>% filter(chr == "11_TNF_PCR_Product") %>% mutate(loci = "TNF")
DC329 <- flip_DNA(DC329)

DC329_ctrl <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 2478) %>% filter(chr == "chr6_(NC_000006)_TNF_long_extraction") %>% mutate(loci = "TNF") %>% 
  filter(condition %in% c("DC329", "DC329-DC324", "DC329-DC332", "DC329-GW134")) %>% 
  mutate(chr = "11_TNF_PCR_Product")
DC329_ctrl <- flip_DNA(DC329_ctrl) 


DC330 <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 215) %>% filter(chr == "12_IL6_PCR_Product") %>% mutate(loci = "IL6")
DC330 <- flip_DNA(DC330)
DC330_ctrl <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 1534) %>% filter(chr == "NC_000007_-_IL6_gene") %>% mutate(loci = "IL6") %>% 
  filter(condition %in% c("DC330", "DC330-DC324", "DC330-DC332", "DC330-GW134"))%>% 
  mutate(chr = "12_IL6_PCR_Product")
DC330_ctrl <- flip_DNA(DC330_ctrl)


DC572 <- MEK1_pileup %>% filter(grepl("DC572", sample)) %>% 
  mutate(n_base = n_base - 130) %>% filter(chr == "01_DNMT1_PCR_Product") %>% mutate(loci = "DNMT1")

DC573 <- MEK1_pileup %>% filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 119) %>% filter(chr == "02_HEK3_PCR_Product") %>% mutate(loci = "HEK3")

DC576 <- MEK1_pileup %>% filter(grepl("DC576", sample)) %>% 
  mutate(n_base = n_base - 251) %>% filter(chr == "05_VEGFA_PCR_Product") %>% mutate(loci = "VEGFA")
DC576 <- flip_DNA(DC576)

DC577 <- MEK1_pileup %>% filter(grepl("DC577", sample)) %>% 
  mutate(n_base = n_base - 263) %>% filter(chr == "06_CD209_PCR_Product") %>% mutate(loci = "CD209")
DC577 <- flip_DNA(DC577)

# Need to change names.
reconcat_all <- rbind(RUNX1, DC203,DC203_ctrl, DC329, DC329_ctrl, DC330, DC330_ctrl, DC572, DC573, DC576, DC577) %>% 
  filter(n_base > -150 & n_base < 150)%>%
  filter(!grepl("GW109", condition))

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
  # Other modes. 
  mutate(GtoC = C/(C+G)*100) %>%
  mutate(GtoT = `T`/(`T`+G)*100) %>%
  mutate(CtoA = A/(C + A)*100) %>% 
  mutate(CtoG = G/(C + G)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) %>% 
  mutate(GtoC = ifelse(ref_base == "G", GtoC, NA)) %>% 
  mutate(GtoT = ifelse(ref_base == "G", GtoT, NA)) %>% 
  mutate(CtoA = ifelse(ref_base == "C", CtoA, NA)) %>% 
  mutate(CtoG = ifelse(ref_base == "C", CtoG, NA)) %>% 
  # Add the control samples. 
  mutate(isControl = ifelse(is.na(helicaseDisplay), 1, 0)) %>% 
  mutate(isBaseline = ifelse(is.na(helicaseDisplay) & is.na(Cas9Display), 1, 0))
  
control_rates <- with_mutation_rate %>% filter(isBaseline  ==  1) %>% 
  group_by(chr, n_base)%>%
  summarise(meanCtoT_ctrl = mean(CtoT), 
            meanGtoA_ctrl = mean(GtoA),
            meanCtoA_ctrl = mean(CtoA), 
            meanCtoG_ctrl = mean(CtoG),
            meanGtoC_ctrl = mean(GtoC), 
            meanGtoT_ctrl = mean(GtoT))

with_mut_and_control <- merge(with_mutation_rate, control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)

# with_mut_background_subtracted <- with_mut_and_control %>% 
#   mutate(CtoT_bg_subtracted = CtoT - meanCtoT_ctrl) %>% 
#   mutate(GtoA_bg_subtracted = GtoA - meanGtoA_ctrl) %>%
#   mutate(CtoA = A/(C + A)) %>%
#   mutate(TtoA = A/(A + `T`))

# ggplot(with_mut_and_control %>% filter(ref_base == "G"), aes(n_base, GtoA, color = helicaseDisplay)) + ylim(0,2) + 
#   geom_point() + facet_grid(helicaseDisplay+Cas9Display ~ loci)

# For each loci, we want to quantify the editing rate.
# allC_bigger_than_10bp <- with_mut_background_subtracted %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
#   filter(!is.na(CtoT_bg_subtracted)) %>%
#   group_by(sample,condition, isControl, loci, helicaseDisplay, Cas9Display) %>% 
#   summarise(meanCtoT_bigger_than_10bp = mean(CtoT_bg_subtracted))
# 
# allC_10bp <- with_mut_background_subtracted %>% filter(ref_base == "C") %>% filter(n_base >= -3 & n_base <= 17) %>% 
#   filter(!is.na(CtoT_bg_subtracted)) %>%
#   group_by(sample) %>%
#   summarise(meanCtoT_spacer = mean(CtoT_bg_subtracted))
# 
# allC_smaller_than_minus10bp <- with_mut_background_subtracted %>% filter(ref_base == "C") %>% filter(n_base < -3 & n_base >=-1000) %>% 
#   filter(!is.na(CtoT_bg_subtracted)) %>%
#   group_by(sample) %>% 
#   summarise(meanCtoT_smaller_than_10bp = mean(CtoT_bg_subtracted))
# 
# allG_10bp <- with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(n_base >= -3 & n_base <= 17) %>% 
#   filter(!is.na(GtoA_bg_subtracted)) %>%
#   group_by(sample) %>% 
#   summarise(meanGtoA_spacer = mean(GtoA_bg_subtracted))
# 
# allG_bigger_than_10bp <- with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(n_base > 17 & n_base <= 1000) %>% 
#   filter(!is.na(GtoA_bg_subtracted)) %>%
#   group_by(sample) %>% 
#   summarise(meanGtoA_bigger_than_10bp = mean(GtoA_bg_subtracted))
# 
# allG_smaller_than_minus10bp <- with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(n_base < -3 & n_base >=-1000) %>% 
#   filter(!is.na(GtoA_bg_subtracted)) %>%
#   group_by(sample) %>% 
#   summarise(meanGtoA_smaller_than_10bp = mean(GtoA_bg_subtracted))

# Not background subtracted. But we still want to filter out the bases that are too high (i.e. real variants). 
with_mutation_rate <- with_mut_and_control %>% 
  filter(meanCtoT_ctrl < 5 | is.na(meanCtoT_ctrl)) %>%
  filter(meanGtoA_ctrl < 5 | is.na(meanGtoA_ctrl))

allC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample,condition, isControl, loci, helicaseDisplay, Cas9Display) %>% 
  summarise(meanCtoT_bigger_than_10bp = mean(CtoT))

allC_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample) %>%
  summarise(meanCtoT_spacer = mean(CtoT))

allC_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample) %>% 
  summarise(meanCtoT_smaller_than_10bp = mean(CtoT))

allG_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  summarise(meanGtoA_spacer = mean(GtoA))

allG_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  summarise(meanGtoA_bigger_than_10bp = mean(GtoA))

allG_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(GtoA)) %>%
  group_by(sample) %>% 
  summarise(meanGtoA_smaller_than_10bp = mean(GtoA))
# GtoC
allGtoC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(GtoC)) %>%
  group_by(sample) %>% 
  summarise(meanGtoC_bigger_than_10bp = mean(GtoC))
allGtoC_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(GtoC)) %>%
  group_by(sample) %>% 
  summarise(meanGtoC_smaller_than_10bp = mean(GtoC))
# GtoT
allGtoT_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(GtoT)) %>%
  group_by(sample) %>% 
  summarise(meanGtoT_bigger_than_10bp = mean(GtoT))
allGtoT_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "G") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(GtoT)) %>%
  group_by(sample) %>% 
  summarise(meanGtoT_smaller_than_10bp = mean(GtoT))
# CtoG
allCtoG_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoG)) %>%
  group_by(sample) %>% 
  summarise(meanCtoG_bigger_than_10bp = mean(CtoG))
allCtoG_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(CtoG)) %>%
  group_by(sample) %>% 
  summarise(meanCtoG_smaller_than_10bp = mean(CtoG))
# CtoA
allCtoA_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoA)) %>%
  group_by(sample) %>% 
  summarise(meanCtoA_bigger_than_10bp = mean(CtoA))
allCtoA_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(CtoA)) %>%
  group_by(sample) %>% 
  summarise(meanCtoA_smaller_than_10bp = mean(CtoA))

# Merge all the stats
merged_stats <- list(allC_10bp, allC_bigger_than_10bp, allC_smaller_than_minus10bp, allG_10bp, allG_bigger_than_10bp, allG_smaller_than_minus10bp,
                     allGtoC_bigger_than_10bp, allGtoC_smaller_than_minus10bp,
                     allGtoT_bigger_than_10bp, allGtoT_smaller_than_minus10bp,
                     allCtoG_bigger_than_10bp, allCtoG_smaller_than_minus10bp,
                     allCtoA_bigger_than_10bp, allCtoA_smaller_than_minus10bp) %>% 
  purrr::reduce(full_join, by = "sample")

merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
  mutate(baseMode = str_extract(mode, "GtoA|CtoT|CtoG|CtoA|GtoT|GtoC"))


# Take average across replicates.
avg_across_reps <- merged_stats_pivot_long %>% 
  filter(!is.na(editRate)) %>%
  group_by(condition, isControl, loci, helicaseDisplay, Cas9Display, baseMode, mode) %>% 
  summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  mutate(baseMode = factor(baseMode, levels = c("GtoA", "CtoT", "CtoG", "CtoA", "GtoT", "GtoC")))

  
##### PLOT NOW. #####
plotting_subset <- avg_across_reps %>%
  filter(baseMode %in% c("GtoA", "CtoT")) %>% 
  filter(helicaseDisplay == "PcrA M6") %>%
  filter( Cas9Display %in% c("nCas9 D10A")| is.na(Cas9Display)) %>% 
  filter(mode != "spacer") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(!grepl("GW109", condition)) %>% 
  mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
  mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer"))) %>% 
  group_by(condition, loci, helicaseDisplay, Cas9Display, modeDisplay) %>% 
  summarise(meanEditRate = sum(meanEditRate), n=n(), sdEditRate = sqrt(sum(sdEditRate^2)))
plotting_subset$Cas9Display[which(!is.na(plotting_subset$Cas9Display))] <- "+nCas9 D10A"
plotting_subset$Cas9Display[which(is.na(plotting_subset$Cas9Display))] <- "-nCas9 D10A"


plotting_subset_points <- merged_stats_pivot_long %>% 
  filter(baseMode %in% c("GtoA", "CtoT")) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!is.na(editRate)) %>%
  filter(helicaseDisplay == "PcrA M6") %>%
  filter( Cas9Display %in% c("nCas9 D10A")| is.na(Cas9Display)) %>% 
  filter(mode != "spacer") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(!grepl("GW109", condition)) %>% 
  mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
  mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer"))) %>% 
  mutate(baseMode = factor(baseMode, levels = c("GtoA", "CtoT"))) %>%
  group_by(sample, condition, loci, helicaseDisplay, Cas9Display, modeDisplay) %>% 
  summarise(editRate = sum(editRate), n = n()) %>% arrange(loci, modeDisplay, Cas9Display)
plotting_subset_points$Cas9Display[which(!is.na(plotting_subset_points$Cas9Display))] <- "+nCas9 D10A"
plotting_subset_points$Cas9Display[which(is.na(plotting_subset_points$Cas9Display))] <- "-nCas9 D10A"
write.csv(plotting_subset_points, file.path(out_dir, "data_fig1g1h.csv"))

# Also look at by mode
plotting_subset_points_by_mode <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("mean\\Sto\\S_", "", mode)) %>% 
  filter(!is.na(editRate)) %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  # filter( Cas9Display %in% c("nCas9 D10A")| is.na(Cas9Display)) %>% 
  # filter(mode != "spacer") %>% 
  # filter(!is.na(helicaseDisplay)) %>% 
  filter(!grepl("GW109", condition)) %>% 
  mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
  mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer"))) %>% 
  mutate(baseMode = factor(baseMode, levels =  c("GtoA", "CtoT", "CtoG", "CtoA", "GtoT", "GtoC"))) %>%
  arrange(loci, modeDisplay, Cas9Display)
# plotting_subset_points_by_mode$Cas9Display[which(!is.na(plotting_subset_points_by_mode$Cas9Display))] <- "+nCas9 D10A"
# plotting_subset_points_by_mode$Cas9Display[which(is.na(plotting_subset_points_by_mode$Cas9Display))] <- "-nCas9 D10A"
write.csv(plotting_subset_points_by_mode, file.path(out_dir, "data_fig1g1h_by_mode.csv"))



g1 <- ggplot(plotting_subset, aes(loci, meanEditRate*10, fill = Cas9Display)) + 
  facet_wrap(~modeDisplay, scales = "free_y") +
  geom_errorbar(aes(ymin=meanEditRate*10-sdEditRate*10/2, ymax=meanEditRate*10+sdEditRate*10), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_point(data=plotting_subset_points, aes(loci, editRate*10, group=Cas9Display), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_startrek() + 
  scale_fill_manual(values = c("#989898", "#cc0b01")) + 
  xlab("Loci") + 
  ylab("Mutations/1000bp") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.9, 0.85)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylim(-0.5,5)
g1

ggsave(file.path(out_dir, "fig1b_range_quantification_updated.pdf"), g1, units = "cm", width = 16, height = 8, dpi = 300)


# ggplot(with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(loci == "VEGFA"), aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) + 
#   # ylim(0,1) + 
#   geom_point() + facet_wrap(~sample, ncol = 6)


## The other plot is going to be that there's no editing if there's only nCas9 or only helicase. 
plotting_subset <- avg_across_reps %>%
  filter(loci == "HEK3") %>% 
  filter(helicaseDisplay == "PcrA M6" | is.na(helicaseDisplay)) %>% 
  filter(Cas9Display %in% c("nCas9 D10A")| is.na(Cas9Display)) %>% 
  filter(mode == "bigger_than_10bp") 


plotting_subset_points <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!is.na(editRate)) %>%
  filter(loci == "HEK3") %>% 
  filter(helicaseDisplay == "PcrA M6" | is.na(helicaseDisplay)) %>% 
  filter(Cas9Display %in% c("nCas9 D10A")| is.na(Cas9Display)) %>% 
  filter(mode == "bigger_than_10bp") %>%
  mutate(baseMode = factor(baseMode, levels = c("GtoA", "CtoT")))


g2<-ggplot(plotting_subset, aes(condition, meanEditRate*10, fill = baseMode)) + 
  geom_errorbar(aes(ymin=meanEditRate*10-sdEditRate*10, ymax=meanEditRate*10+sdEditRate*10), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_point(data=plotting_subset_points, aes(condition, editRate*10, group=baseMode), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  # geom_point(data=plotting_subset_points, aes(loci, editRate, group=baseMode), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_startrek() + 
  xlab("Condition") + 
  ylab("Mutations/1000bp") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.2, 0.8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig1c_dependent_on_Cas9.pdf"), g2, units = "cm", width = 8, height = 8, dpi = 300)


ggplot(plotting_subset, aes(condition, meanEditRate*10, fill = baseMode)) + 
  geom_col(position = position_dodge(0.9, preserve = "single"), color = "black", width = 0.9) + 
  geom_beeswarm(data = plotting_subset_points, aes(x = condition, editRate*10, group = baseMode),
                   dodge.width = 0.9, color = "black", fill = "white", shape=21, cex = 0.5, size = 2.5, method = "swarm") + 
  geom_errorbar(aes(ymin=meanEditRate*10-sdEditRate*10, ymax=meanEditRate*10+sdEditRate*10), width=.2, 
                position=position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = salmon_on_ice) + 
  xlab("Condition") + 
  ylab("Mutations/1000bp")

write.csv(plotting_subset_points, file.path(out_dir, "data_fig1d.csv"))
