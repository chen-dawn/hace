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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/pileup/"

##### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=_\\d+_S\\d+)")) %>%
  select(-filename) %>%
  mutate(base_sum = (A + C + `T` + G )) %>% 
  filter(base_sum > 100)

fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/merged_all_timecourse.csv")

MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/merged_all_timecourse.csv")

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
  mutate(n_base = n_base - 102) %>% filter(chr == "MEK1i1_300bp") %>% mutate(loci = "MAP2K1")

DC573 <- MEK1_pileup %>% # filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 119) %>% filter(chr == "02_HEK3_PCR_Product") %>% mutate(loci = "HEK3")

DC576 <- MEK1_pileup %>% # filter(grepl("DC576", sample)) %>% 
  mutate(n_base = n_base - 251) %>% filter(chr == "05_VEGFA_PCR_Product") %>% mutate(loci = "VEGFA")
DC576 <- flip_DNA(DC576)

DC577 <- MEK1_pileup %>%# filter(grepl("DC577", sample)) %>% 
  mutate(n_base = n_base - 263) %>% filter(chr == "06_CD209_PCR_Product") %>% mutate(loci = "CD209")
DC577 <- flip_DNA(DC577)

# Need to change names.
reconcat_all <- rbind(DC203, DC573, DC576, DC577) %>% 
  filter(n_base > -150 & n_base < 150)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all %>% 
  mutate(helicase = str_extract(condition, "BLM|GSPcrA")) %>% 
  mutate(time = str_extract(condition, "\\d+(?=H)")) %>% 
  mutate(time = ifelse(is.na(time), 0, time)) %>% 
  mutate(time = as.integer(time))


with_mutation_rate <- with_plasmid_annotation %>% 
  filter((A+C+G+`T`)>10000)%>%
  mutate(GtoA = A/(A+G)*100) %>% 
  mutate(CtoT = `T`/(C + `T`)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) 

# ggplot(with_mut_and_control %>% filter(ref_base == "G"), aes(n_base, GtoA, color = helicaseDisplay)) + ylim(0,2) + 
#   geom_point() + facet_grid(helicaseDisplay+Cas9Display ~ loci)

# For each loci, we want to quantify the editing rate.
allC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample, condition, time, loci, helicase) %>% 
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

# sample_exclude <- c("DC576-GW117-Cas9-F05_S107", "DC573-GW117-Cas9-D05_S58", "DC573-GW117-Cas9-E05_S67", "DC573-GW117-Cas9-F05_S76", 
#                     "DC577-GW117-Cas9-A05_S85", "DC577-GW117-Cas9-B05_S94", "DC577-GW117-Cas9-C05_S103",
#                     "DC577-GW119-Cas9-A07_S87", "DC577-GW119-Cas9-B07_S96", "DC577-GW119-Cas9-C07_S105",
#                     "DC573-GW119-Cas9-D07_S60", "DC573-GW119-Cas9-E07_S69", "DC573-GW119-Cas9-F07_S78",
#                     "DC573-pBO111-DC332-D08_S61", "DC573-pBO111-DC332-D09_S62", "DC573-pBO111-DC332-E08_S70",
#                     "DC573-pBO111-DC332-E09_S48", "DC573-pBO111-DC332-E09_S71", "DC573-pBO111-DC332-F08_S56",
#                     "DC573-pBO111-DC332-F08_S79", "DC573-pBO111-DC332-F09_S57", "DC573-pBO111-DC332-F09_S80",
#                     "DC573-GW108-Cas9-D01_S54", "DC573-GW108-Cas9-E01_S63", "DC573-GW108-Cas9-F01_S72")
# 
# exclude_GW110 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW110-Cas9" & mode == "meanGtoA_bigger_than_10bp")
# exclude_GW112 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW112-Cas9" & mode == "meanGtoA_bigger_than_10bp")
# exclude_GW116 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW116-Cas9" & mode == "meanGtoA_bigger_than_10bp")
# exclude_GW118 <- merged_stats_pivot_long %>% filter(condition == "DC573-GW118-Cas9" & mode == "meanGtoA_bigger_than_10bp")
# 
# sample_exclude <- c(sample_exclude, exclude_GW110$sample[c(2,4,6)], exclude_GW112$sample[c(2,4,6)], exclude_GW116$sample[c(2,4,6)], exclude_GW118$sample[c(2,4,6)])
# 
#   
merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
  mutate(baseMode = str_extract(mode, "GtoA|CtoT")) %>% 
  filter(!grepl("GW109", condition)) %>% 
  filter(!grepl("GW111", condition)) %>% 
  filter(!grepl("GW113", condition)) %>% 
  filter(time < 100) %>% 
  mutate(condition_loci = str_extract(condition, "CD209|HEK3")) %>% 
  filter(!sample %in% c("HEK3_GSPcrA_96H_1_S285", "HEK3_GSPcrA_96H_2_S286", "HEK3_GSPcrA_96H_3_S287"))


# Take average across replicates.
avg_across_reps <- merged_stats_pivot_long %>% 
  filter(!is.na(editRate)) %>%
  group_by(loci, helicase, time, baseMode, mode) %>% 
  summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) 

plotting_subset <- merged_stats_pivot_long %>% 
  filter(grepl("bigger_than_10bp", mode)) %>% 
  filter(condition_loci == "HEK3") %>% 
  filter(loci == "HEK3") %>% 
  filter(helicase == "GSPcrA" | is.na(helicase)) %>% 
  group_by(sample, condition, loci, helicase, time) %>% 
  summarise(editRateSum = sum(editRate)) %>% 
  group_by(condition, loci, helicase, time) %>% 
  summarise(meanEditRate = mean(editRateSum*10), sdEditRate = sd(editRateSum*10)) %>%
  mutate(time_string = as.character(time)) 

f2 <- ggplot(plotting_subset, aes(time_string, meanEditRate, group = 1)) + 
  geom_line(color ="#2185C5") + 
  geom_errorbar(aes(ymin=meanEditRate-sdEditRate, ymax=meanEditRate+sdEditRate), color = "#2185C5", width = 0.1) + 
  geom_point(color ="#2185C5") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("Time (h)") + 
  ylab("Mutations/1000bp") 
ggsave(file.path(out_dir, "fig1g_supp_mut_per_1000bp.pdf"), f2, units = "cm", width = 11, height = 8, dpi = 300)


# plotting_subset_points <- merged_stats_pivot_long %>%  
#   mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
#   mutate(mode = gsub("meanGtoA_", "", mode)) %>%
#   filter(mode == "bigger_than_10bp") %>%
#   # filter(baseMode == "GtoA") %>%
#   filter(loci == "HEK3") %>%
#   filter(helicase == "GSPcrA") 
# 
# ggplot(plotting_subset_points, aes(time, editRate, color = baseMode)) + geom_point() + facet_wrap(~condition_loci)
# 
# group_by(condition, time) %>%
#   summarise(mean_mut_per_read = mean(mut_per_read), sd_mut_per_read = sd(mut_per_read, na.rm = T)) %>%
#   filter(grepl("GSPcrA", condition)| grepl("WT", condition)) %>%
#   mutate(time_string = as.character(time)) 
# 


# 
# 
# 
# ##### PLOT NOW. #####
# plotting_subset <- avg_across_reps %>%
#   filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
#   filter(UGIDisplay == "Yes") %>%
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter(Cas9Display %in% c("nCas9 D10A")) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay)) %>% 
#   filter(baseMode == "GtoA")
#   # filter(!grepl("GW109", condition)) %>% 
#   #mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
#   # mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))
# 
# plotting_subset_points <- merged_stats_pivot_long %>% 
#   mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
#   mutate(mode = gsub("meanGtoA_", "", mode)) %>%
#   filter(!grepl("GW109", condition)) %>% 
#   filter(!is.na(editRate)) %>%
#   filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
#   filter(UGIDisplay == "Yes") %>%
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter( Cas9Display %in% c("nCas9 D10A")) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay)) %>% 
#   filter(baseMode == "GtoA")
#   
# g1 <- ggplot(plotting_subset, aes(loci, meanEditRate, fill = helicaseDisplay, group = helicaseDisplay)) + 
#   # facet_grid(Cas9Display~modeDisplay) +
#   geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   geom_point(data=plotting_subset_points, aes(loci, editRate, group=helicaseDisplay), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
#   theme_bw()+ 
#   scale_fill_startrek() + 
#   xlab("Loci") + 
#   ylab("Mean Edits Per kb") + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = c(0.9, 0.8)) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file.path(out_dir, "fig2a_different_helicases_V2.pdf"), g1, units = "cm", width = 14, height = 10, dpi = 300)
# 
# 
# # ggplot(with_mut_background_subtracted %>% filter(ref_base == "G") %>% filter(loci == "VEGFA"), aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) + 
# #   # ylim(0,1) + 
# #   geom_point() + facet_wrap(~sample, ncol = 6)
# 
# 
# ## The other plot is going to be that UGI is important. 
# plotting_subset <- avg_across_reps %>%
#   filter(loci %in% c("HEK3")) %>% 
#   # filter(UGIDisplay == "Yes") %>%
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter(Cas9Display %in% c("nCas9 D10A")) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay)) %>% 
#   filter(baseMode == "GtoA") %>% 
#   mutate(UGIDisplay = factor(UGIDisplay, levels = c("Yes", "No"))) 
# # filter(!grepl("GW109", condition)) %>% 
# #mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# # mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))
# 
# plotting_subset_points <- merged_stats_pivot_long %>% 
#   mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
#   mutate(mode = gsub("meanGtoA_", "", mode)) %>%
#   filter(!is.na(editRate)) %>%
#   filter(loci %in% c("HEK3")) %>% 
#   # filter(UGIDisplay == "Yes") %>%
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter( Cas9Display %in% c("nCas9 D10A")) %>% 
#   filter(!grepl("GW109", condition)) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay)) %>% 
#   filter(baseMode == "GtoA") %>% 
#   mutate(UGIDisplay = factor(UGIDisplay, levels = c("Yes", "No")))
# 
# g2 <- ggplot(plotting_subset, aes(helicaseDisplay, meanEditRate, fill = UGIDisplay, group = interaction(UGIDisplay, loci))) + 
#   # facet_grid(Cas9Display~modeDisplay) +
#   geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   geom_point(data=plotting_subset_points, aes(helicaseDisplay, editRate, group=interaction(UGIDisplay, loci)), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
#   theme_bw()+ 
#   scale_fill_startrek() + 
#   xlab("Helicase") + 
#   ylab("Mean G>A Edit Rate") + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = c(0.9, 0.85)) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file.path(out_dir, "fig2b_with_without_UGI.pdf"), g2, units = "cm", width = 10, height = 10, dpi = 300)
# 
# 
# # This plot is going to be different Cas9 how it affects it. 
# plotting_subset <- avg_across_reps %>%
#   filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
#   filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>% 
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay))
#   # filter(baseMode == "GtoA")
# # filter(!grepl("GW109", condition)) %>% 
# #mutate(modeDisplay = ifelse(mode == "smaller_than_10bp", "Before sgRNA Spacer", "After sgRNA Spacer")) %>% 
# # mutate(modeDisplay = factor(modeDisplay, levels = c("Before sgRNA Spacer", "After sgRNA Spacer")))
# 
# plotting_subset_points <- merged_stats_pivot_long %>% 
#   mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
#   mutate(mode = gsub("meanGtoA_", "", mode)) %>%
#   filter(!grepl("GW109", condition)) %>% 
#   filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>% 
#   filter(!is.na(editRate)) %>%
#   filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3")) %>% 
#   filter(UGIDisplay == "Yes") %>%
#   # filter(helicaseDisplay == "PcrA M6") %>%
#   filter(helicaseDisplay != "RepX") %>%
#   filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
#   filter(mode == "bigger_than_10bp") %>% 
#   filter(!is.na(helicaseDisplay))  
#   # filter(baseMode == "GtoA")
# 
# g3 <- ggplot(plotting_subset, aes(helicaseDisplay, meanEditRate, fill = loci, group = loci)) + 
#   # facet_grid(Cas9Display~modeDisplay) +
#   geom_errorbar(aes(ymin=meanEditRate-sdEditRate/2, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   geom_point(data=plotting_subset_points, aes(helicaseDisplay, editRate, group=loci), position = position_dodge(width = 0.9), size=1, show.legend=FALSE) + 
#   theme_bw()+ 
#   facet_grid(baseMode~Cas9Display) + 
#   scale_fill_npg() + 
#   xlab("Helicase") + 
#   ylab("Mean G>A Edit Rate") + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = c(0.1, 0.82)) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file.path(out_dir, "fig2c_modeCas9.pdf"), g3, units = "cm", width = 16, height = 14, dpi = 300)
# 
# ##### Heatmap #####
# avg_across_reps_GtoA <- avg_across_reps %>% # filter(baseMode == "GtoA") %>% 
#   filter(helicaseDisplay %in% c("PcrA M6", "BLM", "Ns3h", "PcrA")) %>%
#   filter(UGIDisplay == "Yes") %>%
#   filter(Cas9Display %in% c("nCas9 D10A", "nCas9 H840A")) %>% 
#   mutate(Cas9Display = factor(Cas9Display, levels = c("nCas9 D10A", "nCas9 H840A"))) %>%
#   mutate(mode == "bigger_than_GtoA")
# 
# g4 <- ggplot(avg_across_reps_GtoA, aes(loci, Cas9Display, fill = meanEditRate)) + 
#   geom_tile(color = "white") + facet_grid(helicaseDisplay~baseMode) + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   scale_fill_material("blue") + coord_equal() + 
#   theme(panel.spacing = unit(0, "lines")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file.path(out_dir, "fig2f_heatmap.pdf"), g4, units = "cm", width = 18, height = 20, dpi = 300)
# 
# 
# 
# 
# 
# ##### Look at number of mismatches different times #####
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/NM_tag/"
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="NM_counter.tsv", full.names = T)
# 
# all_files <- vroom(MEK1_pileup_filenames, id = "filename", col_names = F)
# all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
#   mutate(condition = str_extract(sample, "^.+(?=_\\d+_S\\d+)")) %>%
#   select(-filename) 
# 
# names(all_files_df) <- c("num_mismatch", "count", "sample", "condition")
# fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/NM_tag/merged_all_NM_tag.csv")
# 
# NM_tag <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/NM_tag/merged_all_NM_tag.csv")
# NM_tag <- NM_tag %>% mutate(num_mismatch = ifelse(num_mismatch > 6, 6, num_mismatch)) %>% 
#   group_by(sample, condition, num_mismatch) %>% 
#   summarise(count = sum(count)) %>% 
#   mutate(num_mismatch = as.character(num_mismatch)) %>%
#   mutate(helicase = str_extract(condition, "BLM|GSPcrA")) %>% 
#   mutate(time = str_extract(condition, "\\d+(?=H)")) %>% 
#   mutate(time = as.integer(time)) %>% 
#   mutate(loci = str_extract(condition, "^[^_]*")) %>% 
#   mutate(rep = str_extract(sample, "\\d(?=_S\\d+)"))
# sample_sums <- NM_tag %>% group_by(sample) %>% summarise(total_count = sum(count))
# NM_tag <- merge(NM_tag, sample_sums, by = "sample") %>% as_tibble()
# 
# # If there's more than 10 mismatch, just group them. 
# NM_tag_processed <- NM_tag %>% mutate(normalized_count = count/total_count*100)
# 
# ggplot(NM_tag_processed %>% filter(!is.na(helicase)& total_count > 10000), aes(time, normalized_count, fill = forcats::fct_rev(num_mismatch))) + 
#   geom_bar(stat = "identity") + # scale_fill_brewer() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(loci+helicase~rep)
# 
# 
