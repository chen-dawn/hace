library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)

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
 
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/"
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)
# 
# all_files <- vroom(MEK1_pileup_filenames, id = "filename")
# all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>% 
#   mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>% 
#   select(-filename)
# 
# all_files_df <- all_files_df %>% 
#   mutate(condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition)) %>% 
#   mutate(condition = gsub("_", "-", condition))
# 
# fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")

c <- reconcat_all %>% group_by(sample) %>% summarise(avg = mean(base_sum))

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
with_plasmid_annotation <- with_names_edited %>% 
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
  mutate(isControl = ifelse(condition %in% c("RUNX1-Ctrl", "DC203-Control", "DC329-Ctrl", "DC330-Ctrl", "DC331-Ctrl", "DC573-pBO101", "DC576-pBO101", "DC577-pBO101"), 1, 0))
  
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

ggplot(with_mut_background_subtracted %>% filter(ref_base.x == "G"), aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) + 
  ylim(0,1) + 
 geom_point() + facet_grid(helicaseDisplay+Cas9Display ~ loci) + xlim(-1000,1000) 

ggplot(with_mut_background_subtracted %>% filter(ref_base.x == "G") %>% filter(loci == "VEGFA"), 
       aes(n_base, GtoA_bg_subtracted, color = Cas9Display)) + 
  # ylim(0,10) + 
  geom_point() + facet_wrap(~ sample) + xlim(-1000,1000) 


avg_mut_per_base <- with_mut_background_subtracted %>% 
  pivot_longer(cols = contains("bg_subtracted"), names_to = "editMode") %>% 
  group_by(chr, n_base, ref_base.x, editMode, condition, loci, helicaseDisplay, Cas9Display) %>% 
  summarise(editRate = mean(value), sd = sd(value)) 

avg_mut_per_base_filtered <- avg_mut_per_base %>% 
  filter(ref_base.x %in% c("C", "G"))

##### PLOT NOW. #####
# plotting_subset <- avg_mut_per_base_filtered %>% 
#   filter(loci == "HEK3") 
#   # filter(condition %in% c("DC573-pBO101-DC332"))
# # plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW113-GW133"))
# 
# p1a_range <- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
#   geom_pointrange(fill='#01aeef', color='grey', shape=21, fatten = 0, size = 0.08) +
#   geom_point(size = 1) + 
#   ylim(0,2) +
#   scale_x_continuous(breaks = round(seq(-1600, 1600, by = 400),1), limits = c(-1000, 1000)) + 
#   theme_bw() + 
#   scale_color_npg() + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   # theme(legend.position = "none") +
#   geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
#   ylab("Mutation rate (%)") + xlab("Base Position") + theme(legend.position = c(0.9, 0.87))

plotting_subset <- avg_mut_per_base_filtered %>% 
  filter(loci == "HEK3") %>% 
  filter(!is.na(Cas9Display)) %>% 
  mutate(Cas9Display = factor(Cas9Display, levels = c("nCas9 D10A", "nCas9 H840A", "dCas9")))
# plotting_subset <- avg_mut_per_base_filtered %>% filter(condition %in% c("DC329-GW113-GW133"))
p2d<- ggplot(data = plotting_subset, aes(n_base, editRate, color = editMode, ymin=editRate-sd, ymax=editRate+sd)) + 
  geom_pointrange(fill='#01aeef', color='grey', shape=21, fatten = 0, size = 0.08) +
  geom_point(size = 1) + 
  ylim(0,2) +
  scale_x_continuous(breaks = round(seq(-1600, 1600, by = 400),1), limits = c(-1000, 1000)) + 
  theme_bw() + 
  scale_color_npg() + 
  facet_grid(helicaseDisplay~Cas9Display) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2 , size=0.5) + 
  ylab("Mutation rate (%)") + xlab("Base Position") # + theme(legend.position = c(0.9, 0.87))
  # facet_zoom(ylim = c(0, 0.5), xlim = c(-1000,1000), zoom.size = 0.5, horizontal = FALSE, show.area = F)
# g <- arrangeGrob(p1a_range, p1a_range2, ncol = 1)

ggsave(file.path(out_dir, "fig2d_range_individual.pdf"), p2d, units = "cm", width = 26, height = 20, dpi = 300)


# 
# 
# with_mut_background_subtracted <- with_mut_background_subtracted %>% mutate(ref_base = ref_base.x)
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
# 
# # Merge all the stats
# merged_stats <- list(allC_10bp, allC_bigger_than_10bp, allC_smaller_than_minus10bp, allG_10bp, allG_bigger_than_10bp, allG_smaller_than_minus10bp) %>% 
#   reduce(full_join, by = "sample")
# 
# merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
#   mutate(baseMode = str_extract(mode, "GtoA|CtoT"))
# 
# 
# # Take average across replicates.
# avg_across_reps <- merged_stats_pivot_long %>% 
#   filter(!is.na(editRate)) %>%
#   group_by(condition, isControl, loci, helicaseDisplay, Cas9Display, baseMode, mode) %>% 
#   summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
#   mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
#   mutate(mode = gsub("meanGtoA_", "", mode)) 
# 
# 
# ##### PLOT NOW. #####
# plotting_subset <- avg_across_reps  %>% # filter( Cas9Display == "nCas9 D10A") %>% 
#   filter(mode != "spacer") %>% 
#   filter(!is.na(helicaseDisplay)) %>% 
#   filter(!grepl("GW109", condition))
# ggplot(plotting_subset, aes(loci, meanEditRate, fill = baseMode)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(Cas9Display+ mode~helicaseDisplay) +
#   geom_errorbar(aes(ymin=meanEditRate-sdEditRate, ymax=meanEditRate+sdEditRate), width=.2, position=position_dodge(.9)) +
#   theme_bw()+ scale_fill_npg() + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = "none") 
