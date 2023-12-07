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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt17_SF3B1_PE/"

##### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename)

# fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt17_SF3B1_PE/merged_all_PE.csv")
# 
# MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt17_SF3B1_PE/merged_all_PE.csv")
MEK1_pileup <- all_files_df

# Merge base sums. 
MEK1_pileup <- MEK1_pileup %>% 
  mutate(sample_without_S = gsub("_S\\d+", "", sample)) %>%
  dplyr::group_by(sample, condition, chr, n_base, ref_base) %>%
  dplyr::summarise(A = sum(A), C = sum(C), G = sum(G), `T` = sum(`T`))

MEK1_pileup <- MEK1_pileup %>% 
  mutate(base_sum = (A + C + `T` + G )) %>% 
  filter(base_sum > 3000)

PE_locations <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt17_SF3B1_PE/SF3B1_PE_locations.csv") %>% 
  mutate(condition = paste0("PE-", plasmidNum))

MEK1_pileup <- merge(MEK1_pileup, PE_locations, by = "condition", all.x = T)


with_mut_rate <- MEK1_pileup %>%
  mutate(GtoA = A/(A+G)) %>% 
  mutate(CtoT = `T`/(C + `T`)) %>% 
  mutate(TtoC = C/(`T` + C)) %>% 
  mutate(AtoG = G/(A+G)) %>% 
  mutate(CtoA = A/(C + A)) %>%
  mutate(TtoA = A/(A + `T`))


# Add mut rate. 
high_mut_rate_all <- with_mut_rate %>% mutate(base_sum = (A+C+G+`T`)) %>% 
  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
       C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
       G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
       T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum)) %>%
  mutate(mut_perc= A_perc + C_perc + G_perc + T_perc) %>% 
  filter(mut_perc > 0.003)
ggplot(high_mut_rate_all %>% filter(n_base < 4000), aes(n_base, mut_perc)) + geom_point() + facet_wrap(~PlasmidFullName, scales = 'free_x') + scale_y_log10()



high_GtoA <- with_mut_rate %>% filter(ref_base == "G" & GtoA > 0.005) %>% 
  mutate(AtoG = 0) %>% 
  mutate(CtoT = 0) %>% 
  mutate(TtoC = 0) 
high_AtoG <- with_mut_rate %>% filter(ref_base == "A" & AtoG > 0.005) %>%
  mutate(GtoA = 0) %>% 
  mutate(CtoT = 0) %>% 
  mutate(TtoC = 0) 
high_CtoT <- with_mut_rate %>% filter(ref_base == "C" & CtoT > 0.005) %>%
  mutate(AtoG = 0) %>% 
  mutate(GtoA = 0) %>% 
  mutate(TtoC = 0) 
high_TtoC <- with_mut_rate %>% filter(ref_base == "T" & TtoC > 0.01) %>%
  mutate(AtoG = 0) %>% 
  mutate(CtoT = 0) %>% 
  mutate(GtoA = 0) 

high_mut_only <- rbind(high_GtoA, high_AtoG, high_CtoT, high_TtoC) %>% as_tibble()

high_mut_only$max_mut <- as.numeric(apply(high_mut_only, 1, function(x){return(max(x[c("AtoG", "GtoA", "CtoT", "TtoC")]))}))
ggplot(high_mut_only, aes(n_base, max_mut)) + geom_point() + facet_wrap(~PlasmidFullName, scales = "free_x")+ scale_y_log10()
# read plasmid locations




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
reconcat_all <- rbind(RUNX1, DC572, DC573, DC576, DC577) %>% 
  filter(n_base > -150 & n_base < 150)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
without_UGI = c("GW108", "GW110", "GW112", "GW114", "GW116", "GW118", "pBO111")
with_plasmid_annotation <- reconcat_all %>% 
  mutate(helicase = str_extract(condition, "DC333|DC334|DC335|pBO101|pBO111|GW109|GW111|GW113|GW108|GW110|GW112|GW114|GW115|GW116|GW117|GW118|GW119|pBO127|DC649|DC391|pBO173")) %>% 
  mutate(Cas9 = str_extract(condition, "DC332|DC324|GW134|GW133|DC327|DC326|dCas9|Cas9")) %>% 
  mutate(helicaseDisplay = helicase) %>% 
  mutate(Cas9Display = Cas9) %>% 
  mutate(UGIDisplay = ifelse(helicase %in% without_UGI, "No", "Yes")) %>% 
  mutate(deaminase = str_extract(condition, "pBO101|DC649|DC391|pBO127|pBO173")) %>%
  mutate(deaminaseDisplay = deaminase)


with_plasmid_annotation$deaminaseDisplay <- str_replace(with_plasmid_annotation$deaminaseDisplay, "pBO101", "AID")
with_plasmid_annotation$deaminaseDisplay <- str_replace(with_plasmid_annotation$deaminaseDisplay, "DC649", "tadA")
with_plasmid_annotation$deaminaseDisplay <- str_replace(with_plasmid_annotation$deaminaseDisplay, "DC391", "tadDE")
with_plasmid_annotation$deaminaseDisplay <- str_replace(with_plasmid_annotation$deaminaseDisplay, "pBO127", "rAPOBEC1")
with_plasmid_annotation$deaminaseDisplay <- str_replace(with_plasmid_annotation$deaminaseDisplay, "pBO173", "TadDE2")



with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC333", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW109|GW108", "BLM")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC334", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW111|GW110", "Ns3h")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "DC335|GW113|GW112", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "GW113|GW112", "PcrA")
with_plasmid_annotation$helicaseDisplay <- str_replace(with_plasmid_annotation$helicaseDisplay, "pBO101|pBO111|pBO127|DC649|DC391|pBO173", "PcrA M6")
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
  mutate(TtoC = C/(`T` + C)*100) %>% 
  mutate(AtoG = G/(A+G)*100) %>% 
  mutate(CtoA = A/(C + A) * 100) %>%
  mutate(TtoA = A/(A + `T`)* 100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA)) %>% 
  mutate(AtoG = ifelse(ref_base == "A", AtoG, NA)) %>% 
  mutate(TtoC = ifelse(ref_base == "T", TtoC, NA)) %>% 
  # Add the control samples. 
  mutate(isControl = ifelse(is.na(helicaseDisplay), 1, 0)) %>% 
  mutate(isBaseline = ifelse(is.na(helicaseDisplay) & is.na(Cas9Display), 1, 0))
  
control_rates <- with_mutation_rate %>% filter(isBaseline  ==  1) %>% 
  group_by(chr, n_base)%>%
  summarise(meanCtoT_ctrl = mean(CtoT), 
            meanGtoA_ctrl = mean(GtoA),
            meanTtoC_ctrl = mean(TtoC),
            meanAtoG_ctrl = mean(AtoG))

with_mut_and_control <- merge(with_mutation_rate, control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)

# Remove the background subtraction.
with_mut_background_subtracted <- with_mut_and_control %>% 
  filter(GtoA < 5) %>%
  filter(CtoT < 5) %>%
  filter(AtoG < 5) %>%
  filter(TtoC < 5) 
  
# For each loci, we want to quantify the editing rate.
allC_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "C") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(CtoT)) %>%
  group_by(sample,condition, isControl, loci, helicaseDisplay, Cas9Display, UGIDisplay, deaminaseDisplay) %>% 
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

allT_10bp <- with_mutation_rate %>% filter(ref_base == "T") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(TtoC)) %>%
  group_by(sample) %>% 
  summarise(meanTtoC_spacer = mean(TtoC))

allT_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "T") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(TtoC)) %>%
  group_by(sample) %>% 
  summarise(meanTtoC_bigger_than_10bp = mean(TtoC))

allT_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "T") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(TtoC)) %>%
  group_by(sample) %>% 
  summarise(meanTtoC_smaller_than_10bp = mean(TtoC))

allA_10bp <- with_mutation_rate %>% filter(ref_base == "A") %>% filter(n_base >= -3 & n_base <= 17) %>% 
  filter(!is.na(AtoG)) %>%
  group_by(sample) %>% 
  summarise(meanAtoG_spacer = mean(AtoG))

allA_bigger_than_10bp <- with_mutation_rate %>% filter(ref_base == "A") %>% filter(n_base > 17 & n_base <= 1000) %>% 
  filter(!is.na(AtoG)) %>%
  group_by(sample) %>% 
  summarise(meanAtoG_bigger_than_10bp = mean(AtoG))

allA_smaller_than_minus10bp <- with_mutation_rate %>% filter(ref_base == "A") %>% filter(n_base < -3 & n_base >=-1000) %>% 
  filter(!is.na(AtoG)) %>%
  group_by(sample) %>% 
  summarise(meanAtoG_smaller_than_10bp = mean(AtoG))

# Merge all the stats
merged_stats <- list(allC_10bp, allC_bigger_than_10bp, allC_smaller_than_minus10bp, 
                     allG_10bp, allG_bigger_than_10bp, allG_smaller_than_minus10bp,
                     allA_10bp, allA_bigger_than_10bp, allA_smaller_than_minus10bp,
                     allT_10bp, allT_bigger_than_10bp, allT_smaller_than_minus10bp) %>% 
  purrr::reduce(full_join, by = "sample")
merged_stats_pivot_long <- merged_stats %>% pivot_longer(cols = starts_with("mean"), names_to = "mode", values_to = "editRate") %>% 
  mutate(baseMode = str_extract(mode, "GtoA|CtoT|TtoC|AtoG")) %>% 
  filter(!grepl("GW109", condition)) %>% 
  filter(!grepl("GW111", condition)) %>% 
  filter(!grepl("GW113", condition)) 


# Take average across replicates.
avg_across_reps <- merged_stats_pivot_long %>% 
  filter(!is.na(editRate)) %>%
  group_by(loci, helicaseDisplay, Cas9Display, deaminaseDisplay, baseMode, mode) %>% 
  summarise(meanEditRate = mean(editRate), sdEditRate = sd(editRate)) %>% 
  mutate(mode = gsub("meanCtoT_", "", mode)) %>% 
  mutate(mode = gsub("meanGtoA_", "", mode)) %>%
  mutate(mode = gsub("meanTtoC_", "", mode)) %>%
  mutate(mode = gsub("meanAtoG_", "", mode)) 
  
ggplot(avg_across_reps %>% filter(mode %in% c("bigger_than_10bp")), aes(deaminaseDisplay, meanEditRate)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(baseMode~loci, scales = "free")

write.csv(merged_stats_pivot_long, file.path(out_dir, "data_fig2m.csv"))

##### PLOT NOW. #####
plotting_subset <- avg_across_reps %>%
  filter(loci %in% c("MAP2K1", "RUNX1", "CD209", "VEGFA", "HEK3", "DNMT1")) %>% 
  filter(UGIDisplay == "Yes") %>%
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
write.csv(plotting_subset_points, file.path(out_dir, "data_fig2m.csv"))

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
ggsave(file.path(out_dir, "fig2m_different_deaminase.pdf"), g1, units = "cm", width = 14, height = 10, dpi = 300)

