library(tidyverse)
library(gridExtra)
library("viridis")
library(ggbreak) 
library(ggsci)
library(ggpubr)
library(cowplot)
out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"

get_CtoT_average <- function(MEK1_pileup){
  # Look at C to T editing.
  MEK1_CtoT <- MEK1_pileup %>% 
    filter(ref_base == "C") %>% 
    filter((A+C+G+`T`)>10000)%>%
    mutate(CtoT = `T`/(`T`+C)) %>% 
    mutate(helicase = str_extract(condition, "GW109|GW111|GW113")) %>% 
    mutate(nCas9 = str_extract(condition, "DC324|GW133|DC327")) %>%
    mutate(loci = str_extract(condition, "^[^-]+")) 
  
  MEK1_CtoT_rep_average <- MEK1_CtoT %>% 
    group_by(condition, n_base) %>%
    summarise(meanCtoTperc = mean(CtoT*100), sd = sd(CtoT*100)) %>%
    mutate(helicase = str_extract(condition, "GW109|GW111|GW113")) %>% 
    mutate(nCas9 = str_extract(condition, "DC324|GW133|DC327")) %>%
    mutate(loci = str_extract(condition, "^[^-]+"))
  
  return(MEK1_CtoT_rep_average)
}


get_GtoA_average <- function(MEK1_pileup){
  MEK1_GtoA <- MEK1_pileup %>% 
    filter(ref_base == "G") %>% 
    filter((A+C+G+`T`)>10000)%>%
    mutate(GtoA = A/(A+G)) %>% 
    mutate(helicase = str_extract(condition, "GW109|GW111|GW113")) %>% 
    mutate(nCas9 = str_extract(condition, "DC324|GW133")) %>%
    mutate(loci = str_extract(condition, "^[^-]+")) 
  
  MEK1_GtoA_rep_average <- MEK1_GtoA %>% 
    group_by(condition, n_base) %>%
    summarise(meanGtoAperc = mean(GtoA*100), sd = sd(GtoA*100)) %>%
    mutate(helicase = str_extract(condition, "GW109|GW111|GW113")) %>% 
    mutate(nCas9 = str_extract(condition, "DC324|GW133|DC327")) %>%
    mutate(loci = str_extract(condition, "^[^-]+"))
  return(MEK1_GtoA_rep_average)
}


### Read in dir.
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/"
# 
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)
# 
# MEK1_pileup <- read_tsv(MEK1_pileup_filenames[1]) %>% 
#   mutate(sample = str_extract(basename(MEK1_pileup_filenames[1]), ".+(S\\d+)")) %>% 
#   mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
# 
# for (i in 2:length(MEK1_pileup_filenames)){
#   tmp <- read_tsv(MEK1_pileup_filenames[i]) %>% 
#     mutate(sample = str_extract(basename(MEK1_pileup_filenames[i]), ".+(S\\d+)")) %>% 
#     mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
#   MEK1_pileup <- rbind(MEK1_pileup, tmp)
# }
# 
# write.csv(MEK1_pileup, "~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/merged_all_A375_selection_DNA_pileups.csv")

MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/merged_all_A375_selection_DNA_pileups.csv")

MEK1_high_coverage <- MEK1_pileup %>% 
  filter((A+C+G+`T`)>10000) %>% 
  mutate(base_sum = (A+C+G+`T`)) %>% 
  mutate(A = A/base_sum*10e6+2) %>% 
  mutate(C = C/base_sum*10e6+2) %>% 
  mutate(G = G/base_sum*10e6+2) %>% 
  mutate(`T` = `T`/base_sum*10e6+2) %>%
  mutate(base_sum = (A+C+G+`T`))

ref_all <- MEK1_high_coverage %>% 
  mutate(A_perc = A/base_sum) %>%
  mutate(C_perc = C/base_sum) %>%
  mutate(G_perc = G/base_sum) %>% 
  mutate(T_perc = `T`/base_sum)

# Select only the S and T drug and also the conditions. Drug1 is the GFP control. 
out_dataframe <-  ref_all %>%  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
                                      C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
                                      G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
                                      T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum)) 

mix2_only <- ref_all %>% filter(grepl("mix2", condition))

mix2_no_drugt <- mix2_only %>% filter(condition == "mix2-day0") %>% mutate(drug = "t") %>% mutate(condition = "mix2t-day0")
mix2_no_drugs <- mix2_only %>% filter(condition == "mix2-day0") %>% mutate(drug = "s") %>% mutate(condition = "mix2s-day0")
mix2_drugs <- mix2_only %>% filter(condition!= "mix2-day0") %>% mutate(drug = str_extract(condition, "(?<=mix2)[a-z]+"))
mix2_all <- rbind(mix2_no_drugt, mix2_no_drugs, mix2_drugs) %>% 
  mutate(day = as.integer(str_extract(condition, "(?<=day)\\d+"))) %>% 
  mutate(GtoA = ifelse(ref_base == "G", A/(G+A), NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", `T`/(C+`T`), NA)) %>%
  mutate(AtoG = ifelse(ref_base == "A", G/(G+A), NA)) %>% 
  mutate(TtoC = ifelse(ref_base == "T", C/(C+`T`), NA))
  

with_no_drug_perc <- mix2_all %>% filter(condition == "mix2t-day0") %>% 
  mutate(drug = "t") %>% 
  mutate(GtoA_ref = GtoA) %>%
  mutate(CtoT_ref = CtoT) %>%
  mutate(AtoG_ref = AtoG) %>%
  mutate(TtoC_ref = TtoC) %>% 
  select(n_base, contains("ref"))

with_drug_folds <- merge(mix2_all %>% filter(!grepl("day0", condition)), with_no_drug_perc, by = c("n_base", "ref_base")) %>% as_tibble()
  
with_drug_folds <- with_drug_folds %>% 
  mutate(GtoA_fold = GtoA/GtoA_ref) %>% 
  mutate(CtoT_fold = CtoT/CtoT_ref) %>% 
  mutate(AtoG_fold = AtoG/AtoG_ref) %>% 
  mutate(TtoC_fold = TtoC/TtoC_ref) 

high_folds_bases <- with_drug_folds %>% 
  filter(day == 20) %>% 
  filter(GtoA_fold > 4 | CtoT_fold > 4 | AtoG_fold > 4 | TtoC_fold > 4) %>% 
  select(condition, drug, day, n_base, ref_base, contains("fold")) 

with_drug_folds_subset <- with_drug_folds %>% filter(n_base %in% high_folds_bases$n_base) %>% 
select(condition, drug, day, n_base, ref_base, contains("fold")) %>% 
pivot_longer(cols = contains("fold")) %>% arrange(n_base)

ggplot(with_drug_folds_subset, aes(day, value, group = n_base)) + geom_point() + geom_line() + 
  facet_grid(drug~name) + theme_bw()

base_with_high_fold_change_t <- mix2_all %>% 
  filter(condition %in% c("mix2t-day0", "mix2t-day20")) %>% 
  select(condition, n_base, ref_base, GtoA) %>% 
  pivot_wider(names_from = "condition", values_from = "GtoA") %>% 
  mutate(GtoA_fold = `mix2t-day20`/`mix2t-day0`) %>% 
  filter(ref_base == "G") %>% 
  filter(GtoA_fold > 4)

base_with_high_fold_change_s <- mix2_all %>% 
  filter(condition %in% c("mix2s-day0", "mix2s-day20")) %>% 
  select(condition, n_base, ref_base, GtoA) %>% 
  pivot_wider(names_from = "condition", values_from = "GtoA") %>% 
  mutate(GtoA_fold = `mix2s-day20`/`mix2s-day0`) %>% 
  filter(ref_base == "G") %>% 
  filter(GtoA_fold > 4)

mix2_all_t_subset <- mix2_all %>% filter(drug == "t") %>% 
  filter(n_base %in% base_with_high_fold_change_t$n_base)
mix2_all_s_subset <- mix2_all %>% filter(drug == "s") %>% 
  filter(n_base %in% base_with_high_fold_change_s$n_base)
mix2_all_GtoA_subset <- rbind(mix2_all_t_subset, mix2_all_s_subset)

ggplot(mix2_all_GtoA_subset %>% filter(ref_base == "G"), aes(day, log(GtoA), group = n_base, label = n_base)) + 
  geom_line() + geom_label() + facet_wrap(~drug) + theme_bw() 


remove_ref_base_coverage <- function(x){
  print(x)
  if (x['ref_base'] == "A"){
    x['A'] = 0
  }
  return (x)
}



mix2_all$exon_num <- NA
mix2_all$exon_num[which(mix2_all$n_base>48090 & mix2_all$n_base < 48400)] <- "Exon 2"
mix2_all$exon_num[which(mix2_all$n_base>49810 & mix2_all$n_base < 50100)] <- "Exon 3"
mix2_all$exon_num[which(mix2_all$n_base>94640 & mix2_all$n_base < 94960)] <- "Exon 6"

mix2_all$drug[which(mix2_all$drug=="s")] <- "Selumetinib"
mix2_all$drug[which(mix2_all$drug=="t")] <- "Trametinib"

# Edit the n_base
mix2_exon2 <- mix2_all %>% filter(n_base >= 48116 & n_base < 48325) %>% mutate(n_base = n_base - 48116+111)
mix2_exon3 <- mix2_all %>% filter(n_base >= 49835 & n_base < 49981) %>% mutate(n_base = n_base - 49845+322+10)
mix2_exon5 <- mix2_all %>% filter(n_base >= 94844 & n_base < 94968) %>% mutate(n_base = n_base - 94844+599)
mix2_all_updated <- rbind(mix2_exon2, mix2_exon3, mix2_exon5)

d20_only <- mix2_all_updated %>% 
  filter(condition %in% c("mix2t-day20", "mix2s-day20"))  %>%
  filter(!is.na(exon_num))

# g1 <- 
ggplot(d20_only%>% filter(ref_base == "G"), aes(n_base, A_perc, label = n_base)) + geom_point() +
  facet_grid(drug~exon_num, scales = "free_x") + scale_y_log10() + 
  theme_bw() + # geom_label() + 
  scale_color_startrek() + 
  # scale_fill_manual(values = c("#989898", "#cc0b01")) + 
  xlab("Base Position") + 
  ylab("Allele Frequency (G>A)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank())
# ggsave(file.path(out_dir, "fig4a_A375_drug.pdf"), g1, units = "cm", width = 18, height = 10, dpi = 300)

# Select only the S and T drug and also the conditions. Drug1 is the GFP control. 
out_dataframe <-  ref_all %>%  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
                                      C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
                                      G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
                                      T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum))  
  # mutate(A_perc = ifelse(A_perc < 5e-6, 5e-6, A_perc)) %>%
  # mutate(C_perc = ifelse(C_perc < 5e-6, 5e-6, C_perc)) %>%
  # mutate(G_perc = ifelse(G_perc < 5e-6, 5e-6, G_perc)) %>%
  # mutate(T_perc = ifelse(T_perc < 5e-6, 5e-6, T_perc)) 


# Select only the S and T drug and also the conditions. Drug1 is the GFP control. 
out_dataframe_filtered <- out_dataframe %>% filter(condition %in% c("mix2-day0", "mix2s-day20", "mix2t-day20")) %>% 
  group_by(condition, n_base, ref_base)

# CtoT
ref_percentages <- out_dataframe_filtered %>% filter(condition == "mix2-day0") %>% ungroup() %>% 
  select(n_base, A_perc, C_perc, G_perc, T_perc) %>% 
  mutate(A_ref_perc = A_perc) %>% 
  mutate(C_ref_perc = C_perc) %>% 
  mutate(G_ref_perc = G_perc) %>% 
  mutate(T_ref_perc = T_perc) %>% 
  select(n_base, A_ref_perc, C_ref_perc, G_ref_perc, T_ref_perc) 

CtoT <- out_dataframe_filtered %>% select(sample, n_base, ref_base, T_perc, A_perc, C_perc, G_perc) %>% filter(condition != "mix2-day0")
CtoT <- merge(CtoT, ref_percentages, by = "n_base") %>% as_tibble() %>% 
  mutate(CtoT_Fold = T_perc/T_ref_perc) %>% 
  mutate(GtoA_Fold = A_perc/A_ref_perc) %>% 
  mutate(TtoC_Fold = C_perc/C_ref_perc) %>% 
  mutate(AtoG_Fold = G_perc/G_ref_perc) 
# %>%
#   summarise(CtoT_fold_mean = mean(CtoT_Fold),
#             CtoT_fold_sd = sd(CtoT_Fold),
#             GtoA_fold_mean = mean(GtoA_Fold, na.rm = T),
#             GtoA_fold_sd = sd(GtoA_Fold, na.rm = T),
#             TtoC_fold_mean = mean(TtoC_Fold, na.rm = T),
#             TtoC_fold_sd = sd(TtoC_Fold, na.rm = T),
#             AtoG_fold_mean = mean(AtoG_Fold, na.rm = T),
#             AtoG_fold_sd = sd(AtoG_Fold, na.rm = T))

CtoT$CtoT_Fold[which(CtoT$ref_base == "A")] <- NA
CtoT$CtoT_Fold[which(CtoT$ref_base == "T")] <- NA
CtoT$CtoT_Fold[which(CtoT$ref_base == "G")] <- NA
CtoT$GtoA_Fold[which(CtoT$ref_base == "A")] <- NA
CtoT$GtoA_Fold[which(CtoT$ref_base == "T")] <- NA
CtoT$GtoA_Fold[which(CtoT$ref_base == "C")] <- NA
CtoT$TtoC_Fold[which(CtoT$ref_base == "A")] <- NA
CtoT$TtoC_Fold[which(CtoT$ref_base == "C")] <- NA
CtoT$TtoC_Fold[which(CtoT$ref_base == "G")] <- NA
CtoT$AtoG_Fold[which(CtoT$ref_base == "G")] <- NA
CtoT$AtoG_Fold[which(CtoT$ref_base == "T")] <- NA
CtoT$AtoG_Fold[which(CtoT$ref_base == "C")] <- NA

# This is super janky to merge both.
CtoTonly <- CtoT %>% select(condition, n_base, CtoT_Fold) %>% 
  mutate(type = "C>T")
names(CtoTonly) <- c("condition", "n_base", "fold_change", "type")
GtoAonly <- CtoT %>% select(condition, n_base, GtoA_Fold) %>% 
  mutate(type = "G>A")
names(GtoAonly) <- c("condition", "n_base", "fold_change", "type")
TtoConly <- CtoT %>% select(condition, n_base, TtoC_Fold) %>% 
  mutate(type = "T>C")
names(TtoConly) <- c("condition", "n_base", "fold_change", "type")
AtoGonly <- CtoT %>% select(condition, n_base, AtoG_Fold) %>% 
  mutate(type = "A>G")
names(AtoGonly) <- c("condition", "n_base", "fold_change", "type")
merged_both_modes <- rbind(CtoTonly, GtoAonly,TtoConly,AtoGonly)
merged_both_modes$exon_num <- NA
merged_both_modes$exon_num[which(merged_both_modes$n_base>48090 & merged_both_modes$n_base < 48400)] <- "Exon 2"
merged_both_modes$exon_num[which(merged_both_modes$n_base>49810 & merged_both_modes$n_base < 50100)] <- "Exon 3"
merged_both_modes$exon_num[which(merged_both_modes$n_base>94640 & merged_both_modes$n_base < 94960)] <- "Exon 6"
merged_both_modes$condition[which(merged_both_modes$condition=="mix2s-day20")] <- "Selumetinib"
merged_both_modes$condition[which(merged_both_modes$condition=="mix2t-day20")] <- "Trametinib"
merged_both_modes$fold_change[which(!is.finite(merged_both_modes$fold_change))] <- NA
merged_both_modes <- merged_both_modes %>% filter(!is.na(exon_num)) 

ggplot(merged_both_modes %>% mutate(fold_change = ifelse(fold_change >20, 20, fold_change)), aes(n_base, fold_change, color = type)) + geom_point() +
  scale_color_npg() + 
  facet_grid(condition~exon_num, scales = "free_x") + theme_bw() + ylab("Fold Enrichment (log10 scale)") + 
  ylim(0,20) + 
  xlab("Genomic Position") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) 
# ggsave(file="~/Dropbox (Harvard University)/03TRACE/FIguresDraft/SandTdrugSelection.pdf", units = "cm", width = 20, height = 7, dpi = 300) #saves g

ggplot(CtoT, aes(n_base, log10(GtoA_fold_mean), color = condition)) + geom_point() +
  # xlim(48090,48350)
  xlim(49770,50020)
  # xlim(94770,95000) + ylim(0, 10)

# Drug 4 vs 40 as control:
drug4_40 <- out_dataframe_filtered %>% select(n_base, ref_base, A_perc)%>% pivot_wider(names_from = condition, values_from = A_perc)

ggplot(drug4_40, aes(Drug1, Drug5)) + geom_point() + 
  geom_line() + 
  scale_x_log10() + scale_y_log10() +
  scale_x_break(c(48070,48400), scales=1) + 
  scale_x_break(c(49700,50100), scales = 1)  

# out_dataframe <- out_dataframe %>% mutate(index = seq(1, nrow(out_dataframe)))
ggplot(out_dataframe %>% filter(condition %in% c("Drug4", "Drug40", "Sdrug", "Tdrug")), aes(n_base, A_perc)) + 
  geom_point() + facet_wrap(~condition) + xlim(49700,50100) + ggtitle("Exon 2") # + xlim(49700,50100) # xlim(48070,48400)# + 
ggplot(out_dataframe, aes(n_base, C_perc)) + geom_point() + facet_wrap(~condition) + xlim(49700,50100)
ggplot(out_dataframe, aes(n_base, G_perc)) + geom_point() + facet_wrap(~condition) + xlim(49700,50100)
ggplot(out_dataframe, aes(n_base, T_perc)) + geom_point() + facet_wrap(~condition) + xlim(49700,50100)

out_dataframe_merged <- out_dataframe %>% group_by(condition, n_base, ref_base) %>% summarise(A_perc = mean(A_perc), 
                                                            C_perc = mean(C_perc),
                                                            G_perc = mean(G_perc),
                                                            T_perc = mean(T_perc))

write.csv(out_dataframe_merged, "~/Dropbox (Harvard University)/03TRACE/Helicase/A375_long_drug_selection.csv")

