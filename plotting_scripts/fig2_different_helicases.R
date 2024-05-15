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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt02_Multiple_Helicases/"

##### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename)

fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt02_Multiple_Helicases/merged_all_multiple_helicases.csv")

MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt02_Multiple_Helicases/merged_all_multiple_helicases.csv")
MEK1_pileup2 <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt01_Multiple_locations_all_helicase_Cas9/merged_all_amplicon_pileups.csv")

MEK1_pileup <- rbind(MEK1_pileup, MEK1_pileup2)

# Merge base sums. 
MEK1_pileup <- MEK1_pileup %>% 
  mutate(sample_without_S = gsub("_S\\d+", "", sample)) %>%
  dplyr::group_by(sample, condition, chr, n_base, ref_base) %>%
  dplyr::summarise(A = sum(A), C = sum(C), G = sum(G), `T` = sum(`T`))

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
DC203 <- MEK1_pileup %>% # filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 140) %>% filter(chr == "10_MAP2K1_PCR_Product") %>% mutate(loci = "MAP2K1")
DC203_ctrl <- MEK1_pileup %>% # filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 102) %>% filter(chr == "MEK1i1_300bp") %>% mutate(loci = "MAP2K1") %>% 
  filter(condition %in% c("DC203", "DC203-DC324", "DC203-DC332", "DC203-GW134")) %>% 
  mutate(chr = "10_MAP2K1_PCR_Product")

RUNX1 <- MEK1_pileup %>% # filter(grepl("RUNX1", sample)) %>% 
  mutate(n_base = n_base - 146) %>% filter(chr == "RUNX1_300bp") %>% mutate(loci = "RUNX1")

DC329 <- MEK1_pileup %>% # filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 217) %>% filter(chr == "11_TNF_PCR_Product") %>% mutate(loci = "TNF")
DC329 <- flip_DNA(DC329)

DC329_ctrl <- MEK1_pileup %>% # filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 2478) %>% filter(chr == "chr6_(NC_000006)_TNF_long_extraction") %>% mutate(loci = "TNF") %>% 
  filter(condition %in% c("DC329", "DC329-DC324", "DC329-DC332", "DC329-GW134")) %>% 
  mutate(chr = "11_TNF_PCR_Product")
DC329_ctrl <- flip_DNA(DC329_ctrl) 


DC330 <- MEK1_pileup %>% # filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 215) %>% filter(chr == "12_IL6_PCR_Product") %>% mutate(loci = "IL6")
DC330 <- flip_DNA(DC330)
DC330_ctrl <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 1534) %>% filter(chr == "NC_000007_-_IL6_gene") %>% mutate(loci = "IL6") %>% 
  filter(condition %in% c("DC330", "DC330-DC324", "DC330-DC332", "DC330-GW134"))%>% 
  mutate(chr = "12_IL6_PCR_Product")
DC330_ctrl <- flip_DNA(DC330_ctrl)

DC572 <- MEK1_pileup %>% # filter(grepl("DC572", sample)) %>% 
  mutate(n_base = n_base - 130) %>% filter(chr == "01_DNMT1_PCR_Product") %>% mutate(loci = "DNMT1")

DC573 <- MEK1_pileup %>% # filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 119) %>% filter(chr == "02_HEK3_PCR_Product") %>% mutate(loci = "HEK3")

DC576 <- MEK1_pileup %>% # filter(grepl("DC576", sample)) %>% 
  mutate(n_base = n_base - 251) %>% filter(chr == "05_VEGFA_PCR_Product") %>% mutate(loci = "VEGFA")
DC576 <- flip_DNA(DC576)

DC577 <- MEK1_pileup %>%  filter(grepl("DC577", sample)) %>% 
  mutate(n_base = n_base - 263) %>% filter(chr == "06_CD209_PCR_Product") %>% mutate(loci = "CD209")
DC577 <- flip_DNA(DC577)

# Need to change names.
reconcat_all <- rbind(RUNX1, DC203, DC329, DC330, DC572, DC573, DC576, DC577) %>% 
  filter(n_base > -150 & n_base < 150)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
without_UGI = c("GW108", "GW110", "GW112", "GW114", "GW116", "GW118", "pBO111")
with_plasmid_annotation <- reconcat_all %>% 
  mutate(helicase = str_extract(condition, "DC333|DC334|DC335|pBO101|pBO111|GW109|GW111|GW113|GW108|GW110|GW112|GW114|GW115|GW116|GW117|GW118|GW119")) %>% 
  mutate(Cas9 = str_extract(condition, "DC332|DC324|GW134|GW133|DC327|DC326|dCas9|Cas9")) %>% 
  mutate(helicaseDisplay = helicase) %>% 
  mutate(Cas9Display = Cas9) %>% 
  mutate(UGIDisplay = ifelse(helicase %in% without_UGI, "No", "Yes"))

with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC333", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW109|GW108", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC334", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW111|GW110", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC335|GW113|GW112", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW113|GW112", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "pBO101|pBO111", "PcrA M6")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW114|GW115", "RepX")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW116|GW117", "TraL")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW118|GW119", "urvD")

with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "dCas9", "IGNORE")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "Cas9", "nCas9 D10A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC332", "nCas9 D10A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "GW133", "nCas9 D10A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC324", "nCas9 H840A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC326", "nCas9 H840A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "DC327", "nCas9 H840A")
with_plasmid_annotation$Cas9Display <- str_replace(with_plasmid_annotation$Cas9Display, "GW134|IGNORE", "dCas9")




with_mutation_rate <- with_plasmid_annotation %>% 
  filter((A+C+G+`T`)>10000)%>%
  mutate(GtoA = A/(A+G)*100) %>% 
  mutate(CtoT = `T`/(C + `T`)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) %>% 
  # Add the control samples. 
  mutate(isControl = ifelse(is.na(helicaseDisplay), 1, 0)) %>% 
  mutate(isBaseline = ifelse(is.na(helicaseDisplay) & is.na(Cas9Display), 1, 0))
  
control_rates <- with_mutation_rate %>% filter(isBaseline  ==  1) %>% 
  group_by(chr, n_base)%>%
  summarise(meanCtoT_ctrl = mean(CtoT), 
            meanGtoA_ctrl = mean(GtoA))

with_mut_and_control <- merge(with_mutation_rate, control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)

# Remove the background subtraction.
with_mut_background_subtracted <- with_mut_and_control %>% 
  #   mutate(CtoT_bg_subtracted = CtoT - meanCtoT_ctrl) %>% 
  mutate(CtoT_bg_subtracted = CtoT) %>% 
  mutate(GtoA_bg_subtracted = GtoA) %>%
  mutate(CtoA = A/(C + A)) %>%
  mutate(TtoA = A/(A + `T`))

# ggplot(with_mut_and_control %>% filter(ref_base == "G"), aes(n_base, GtoA, color = helicaseDisplay)) + ylim(0,2) + 
#   geom_point() + facet_grid(helicaseDisplay+Cas9Display ~ loci)

# For each loci, we want to quantify the editing rate.
allC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample,condition, isControl, loci, helicaseDisplay, Cas9Display, UGIDisplay) %>% 
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

# Merge all the stats
merged_stats <- list(allC_10bp, allC_bigger_than_10bp, allC_smaller_than_minus10bp, allG_10bp, allG_bigger_than_10bp, allG_smaller_than_minus10bp) %>% 
  purrr::reduce(full_join, by = "sample")

sample_exclude <- c("DC576-GW117-Cas9-F05_S107", "DC573-GW117-Cas9-D05_S58", "DC573-GW117-Cas9-E05_S67", "DC573-GW117-Cas9-F05_S76", 
                    "DC577-GW117-Cas9-A05_S85", "DC577-GW117-Cas9-B05_S94", "DC577-GW117-Cas9-C05_S103",
                    "DC577-GW119-Cas9-A07_S87", "DC577-GW119-Cas9-B07_S96", "DC577-GW119-Cas9-C07_S105",
                    "DC573-GW119-Cas9-D07_S60", "DC573-GW119-Cas9-E07_S69", "DC573-GW119-Cas9-F07_S78",
                    "DC573-pBO111-DC332-D08_S61", "DC573-pBO111-DC332-D09_S62", "DC573-pBO111-DC332-E08_S70",
                    "DC573-pBO111-DC332-E09_S48", "DC573-pBO111-DC332-E09_S71", "DC573-pBO111-DC332-F08_S56",
                    "DC573-pBO111-DC332-F08_S79", "DC573-pBO111-DC332-F09_S57", "DC573-pBO111-DC332-F09_S80",
                    "DC573-GW108-Cas9-D01_S54", "DC573-GW108-Cas9-E01_S63", "DC573-GW108-Cas9-F01_S72")

exclude_GW110 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW110-Cas9" & mode == "meanGtoA_bigger_than_10bp")
exclude_GW112 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW112-Cas9" & mode == "meanGtoA_bigger_than_10bp")
exclude_GW116 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW116-Cas9" & mode == "meanGtoA_bigger_than_10bp")
exclude_GW118 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW118-Cas9" & mode == "meanGtoA_bigger_than_10bp")

sample_exclude <- c(sample_exclude, exclude_GW110$sample[c(2,4,6)], exclude_GW112$sample[c(2,4,6)], exclude_GW116$sample[c(2,4,6)], exclude_GW118$sample[c(2,4,6)])

  
merged_stats <- merged_stats %>% filter(!(sample %in% sample_exclude))
merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
  mutate(baseMode = str_extract(mode, "GtoA|CtoT")) %>% 
  filter(!grepl("GW109", condition)) %>% 
  filter(!grepl("GW111", condition)) %>% 
  filter(!grepl("GW113", condition)) 


# Take average across replicates.
avg_across_reps <- merged_stats_pivot_long %>% 
  filter(!is.na(editRate)) %>%
  group_by(loci, helicaseDisplay, Cas9Display, UGIDisplay, baseMode, mode) %>% 
  summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) 
  
  
##### PLOT NOW. #####
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
  filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter(Cas9Display %in% c("nCas9 D10A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA")
  # filter(!grepl("GW109", condition)) %>% 
  #mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
  # mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))

plotting_subset_points <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!grepl("GW109", condition)) %>% 
  filter(!is.na(editRate)) %>%
  # filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
  filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  # filter( Cas9Display %in% c("nCas9 H840A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA")
write.csv(plotting_subset_points, file.path(out_dir, "data_fig2a.csv"))

newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")

g1 <- ggplot(plotting_subset, aes(loci, meanEditRate, fill = helicaseDisplay, group = helicaseDisplay)) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), color = "black", size = 0.4) + 
  geom_point(data=plotting_subset_points, aes(loci, editRate, group=helicaseDisplay), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_brewer(direction = -1) + 
  # scale_fill_manual(values = newPalette400_rotate) + 
  xlab("Loci") + 
  ylab("Mean Edits Per kb") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.9, 0.8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig2a_different_helicases_V2.pdf"), g1, units = "cm", width = 14, height = 10, dpi = 300)


# ggplot(with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(loci == "VEGFA"), aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) + 
#   # ylim(0,1) + 
#   geom_point() + facet_wrap(~sample, ncol = 6)


## The other plot is going to be that UGI is important. 
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("HEK3")) %>% 
  # filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter(Cas9Display %in% c("nCas9 D10A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA") %>% 
  mutate(UGIDisplay = factor(UGIDisplay, levels = c("Yes", "No"))) 
# filter(!grepl("GW109", condition)) %>% 
#mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))

plotting_subset_points <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!is.na(editRate)) %>%
  filter(loci %in% c("HEK3")) %>% 
  # filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter( Cas9Display %in% c("nCas9 D10A")) %>% 
  filter(!grepl("GW109", condition)) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA") %>% 
  mutate(UGIDisplay = factor(UGIDisplay, levels = c("Yes", "No")))

g2 <- ggplot(plotting_subset, aes(helicaseDisplay, meanEditRate, fill = UGIDisplay, group = interaction(UGIDisplay, loci))) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), color = "black", size = 0.4) + 
  geom_point(data=plotting_subset_points, aes(helicaseDisplay, editRate, group=interaction(UGIDisplay, loci)), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_brewer(palette = "Greys", direction = -1) + 
  #  scale_fill_manual(values = newPalette400_rotate) + 
  xlab("Helicase") + 
  ylab("Mean G>A Edit Rate") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.9, 0.85)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig2b_with_without_UGI.pdf"), g2, units = "cm", width = 10, height = 7, dpi = 300)
write_csv(plotting_subset_points, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_2a_with_without_UGI.csv")

# This plot is going to be different Cas9 how it affects it. 
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
  filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>% 
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay))
  # filter(baseMode == "GtoA")
# filter(!grepl("GW109", condition)) %>% 
#mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))

plotting_subset_points <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!grepl("GW109", condition)) %>% 
  filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>% 
  filter(!is.na(editRate)) %>%
  filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
  filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay))  
  # filter(baseMode == "GtoA")

g3 <- ggplot(plotting_subset, aes(helicaseDisplay, meanEditRate, fill = loci, group = loci)) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_point(data=plotting_subset_points, aes(helicaseDisplay, editRate, group=loci), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
  theme_bw()+ 
  facet_grid(baseMode~Cas9Display) + 
  scale_fill_npg() + 
  xlab("Helicase") + 
  ylab("Mean G>A Edit Rate") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.1, 0.82)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig2c_modeCas9.pdf"), g3, units = "cm", width = 16, height = 14, dpi = 300)

##### Heatmap #####
avg_across_reps_GtoA <- avg_across_reps %>% # filter(baseMode == "GtoA") %>% 
  filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>%
  filter(UGIDisplay == "Yes") %>%
  filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
  mutate(Cas9Display = factor(Cas9Display, levels = c("nCas9 D10A", "nCas9 H840A"))) %>%
  mutate(mode == "bigger_than_GtoA")

g4 <- ggplot(avg_across_reps_GtoA, aes(loci, Cas9Display, fill = meanEditRate)) + 
  geom_tile(color = "white") + facet_grid(helicaseDisplay~baseMode) + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_fill_material("blue") + coord_equal() + 
  theme(panel.spacing = unit(0, "lines")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(out_dir, "fig2f_heatmap.pdf"), g4, units = "cm", width = 18, height = 20, dpi = 300)



###### REDO 1 ONLY 1 LOCI ######
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("CD209")) %>% 
  filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter(Cas9Display %in% c("nCas9 D10A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA")
# filter(!grepl("GW109", condition)) %>% 
#mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))

plotting_subset_points <- merged_stats_pivot_long %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  filter(!grepl("GW109", condition)) %>% 
  filter(!is.na(editRate)) %>%
  filter(loci %in% c("CD209")) %>% 
  filter(UGIDisplay == "Yes") %>%
  # filter(helicaseDisplay == "PcrA M6") %>%
  filter(helicaseDisplay != "RepX") %>%
  filter( Cas9Display %in% c("nCas9 D10A")) %>% 
  filter(mode == "bigger_than_10bp") %>% 
  filter(!is.na(helicaseDisplay)) %>% 
  filter(baseMode == "GtoA")

ggplot(plotting_subset, aes(helicaseDisplay, meanEditRate*10, group = helicaseDisplay)) + 
  # facet_grid(Cas9Display~modeDisplay) +
  geom_errorbar(aes(ymin=(meanEditRate-sdEditRate/2)*10, ymax=(meanEditRate+sdEditRate)*10), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), color = "black", size = 0.4) + 
  geom_dotplot(data=plotting_subset_points, aes(x = helicaseDisplay, y = editRate*10), stackdir = "center", position = "jitter", size=1, show.legend=FALSE) + 
  theme_bw()+ 
  scale_fill_brewer(direction = -1) + 
  ylab("Mean Edits Per kb") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  # theme(legend.position = c(0.9, 0.8)) + 
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file.path(out_dir, "fig2a_different_helicases_V2.pdf"), g1, units = "cm", width = 14, height = 10, dpi = 300)

