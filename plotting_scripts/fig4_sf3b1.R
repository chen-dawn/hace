library(tidyverse)
library(magrittr)
library(dplyr)
library(gridExtra)
library("viridis")
library(ggbreak) 
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(seqinr)
library(Biostrings)
library("PAIRADISE")
out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"


SF3B1_CDS <- "ATGGCGAAGATCGCCAAGACTCACGAAGATATTGAAGCACAGATTCGAGAAATTCAAGGCAAGAAGGCAGCTCTTGATGAAGCTCAAGGAGTGGGCCTCGATTCTACAGGTTATTATGACCAGGAAATTTATGGTGGAAGTGACAGCAGATTTGCTGGATACGTGACATCAATTGCTGCAACTGAACTTGAAGATGATGACGATGACTATTCATCATCTACGAGTTTGCTTGGTCAGAAGAAGCCAGGATATCATGCCCCTGTGGCATTGCTTAATGATATACCACAGTCAACAGAACAGTATGATCCATTTGCTGAGCACAGACCTCCAAAGATTGCAGACCGGGAAGATGAATACAAAAAGCATAGGCGGACCATGATAATTTCCCCAGAGCGTCTTGATCCTTTTGCAGATGGAGGGAAAACCCCTGATCCTAAAATGAATGCTAGGACTTACATGGATGTAATGCGAGAACAACACTTGACTAAAGAAGAACGAGAAATTAGGCAACAGCTAGCAGAAAAAGCTAAAGCTGGAGAACTAAAAGTCGTCAATGGAGCAGCAGCGTCCCAGCCTCCATCAAAACGAAAACGGCGTTGGGATCAAACAGCTGATCAGACTCCTGGTGCCACTCCCAAAAAACTATCAAGTTGGGATCAGGCAGAGACCCCTGGGCATACTCCTTCCTTAAGATGGGATGAGACACCAGGTCGTGCAAAGGGAAGCGAGACTCCTGGAGCAACCCCAGGCTCAAAAATATGGGATCCTACACCTAGCCACACACCAGCGGGAGCTGCTACTCCTGGACGAGGTGATACACCAGGCCATGCGACACCAGGCCATGGAGGCGCAACTTCCAGTGCTCGTAAAAACAGATGGGATGAAACCCCCAAAACAGAGAGAGATACTCCTGGGCATGGAAGTGGATGGGCTGAGACTCCTCGAACAGATCGAGGTGGAGATTCTATTGGTGAAACACCGACTCCTGGAGCCAGTAAAAGAAAATCACGGTGGGATGAAACACCAGCTAGTCAGATGGGTGGAAGCACTCCAGTTCTGACCCCTGGAAAGACACCAATTGGCACACCAGCCATGAACATGGCTACCCCTACTCCAGGTCACATAATGAGTATGACTCCTGAACAGCTTCAGGCTTGGCGGTGGGAAAGAGAAATTGATGAGAGAAATCGCCCACTTTCTGATGAGGAATTAGATGCTATGTTCCCAGAAGGATATAAGGTACTTCCTCCTCCAGCTGGTTATGTTCCTATTCGAACTCCAGCTCGAAAGCTGACAGCTACTCCAACACCTTTGGGTGGTATGACTGGTTTCCACATGCAAACTGAAGATCGAACTATGAAAAGTGTTAATGACCAGCCATCTGGAAATCTTCCATTTTTAAAACCTGATGATATTCAATACTTTGATAAACTATTGGTTGATGTTGATGAATCAACACTTAGTCCAGAAGAGCAAAAAGAGAGAAAAATAATGAAGTTGCTTTTAAAAATTAAGAATGGAACACCACCAATGAGAAAGGCTGCATTGCGTCAGATTACTGATAAAGCTCGTGAATTTGGAGCTGGTCCTTTGTTTAATCAGATTCTTCCTCTGCTGATGTCTCCTACACTTGAGGATCAAGAGCGTCATTTACTTGTGAAAGTTATTGATAGGATACTGTACAAACTTGATGACTTAGTTCGTCCATATGTGCATAAGATCCTCGTGGTCATTGAACCGCTATTGATTGATGAAGATTACTATGCTAGAGTGGAAGGCCGAGAGATCATTTCTAATTTGGCAAAGGCTGCTGGTCTGGCTACTATGATCTCTACCATGAGACCTGATATAGATAACATGGATGAGTATGTCCGTAACACAACAGCTAGAGCTTTTGCTGTTGTAGCCTCTGCCCTGGGCATTCCTTCTTTATTGCCCTTCTTAAAAGCTGTGTGCAAAAGCAAGAAGTCCTGGCAAGCGAGACACACTGGTATTAAGATTGTACAACAGATAGCTATTCTTATGGGCTGTGCCATCTTGCCACATCTTAGAAGTTTAGTTGAAATCATTGAACATGGTCTTGTGGATGAGCAGCAGAAAGTTCGGACCATCAGTGCTTTGGCCATTGCTGCCTTGGCTGAAGCAGCAACTCCTTATGGTATCGAATCTTTTGATTCTGTGTTAAAGCCTTTATGGAAGGGTATCCGCCAACACAGAGGAAAGGGTTTGGCTGCTTTCTTGAAGGCTATTGGGTATCTTATTCCTCTTATGGATGCAGAATATGCCAACTACTATACTAGAGAAGTGATGTTAATCCTTATTCGAGAATTCCAGTCTCCTGATGAGGAAATGAAAAAAATTGTGCTGAAGGTGGTAAAACAGTGTTGTGGGACAGATGGTGTAGAAGCAAACTACATTAAAACAGAGATTCTTCCTCCCTTTTTTAAACACTTCTGGCAGCACAGGATGGCTTTGGATAGAAGAAATTACCGACAGTTAGTTGATACTACTGTGGAGTTGGCAAACAAAGTAGGTGCAGCAGAAATTATATCCAGGATTGTGGATGATCTGAAAGATGAAGCCGAACAGTACAGAAAAATGGTGATGGAGACAATTGAGAAAATTATGGGTAATTTGGGAGCAGCAGATATTGATCATAAACTTGAAGAACAACTGATTGATGGTATTCTTTATGCTTTCCAAGAACAGACTACAGAGGACTCAGTAATGTTGAACGGCTTTGGCACAGTGGTTAATGCTCTTGGCAAACGAGTCAAACCATACTTGCCTCAGATCTGTGGTACAGTTTTGTGGCGTTTAAATAACAAATCTGCTAAAGTTAGGCAACAGGCAGCTGACTTGATTTCTCGAACTGCTGTTGTCATGAAGACTTGTCAAGAGGAAAAATTGATGGGACACTTGGGTGTTGTATTGTATGAGTATTTGGGTGAAGAGTACCCTGAAGTATTGGGCAGCATTCTTGGAGCACTGAAGGCCATTGTAAATGTCATAGGTATGCATAAGATGACTCCACCAATTAAAGATCTGCTGCCTAGACTCACCCCCATCTTAAAGAACAGACATGAAAAAGTACAAGAGAATTGTATTGATCTTGTTGGTCGTATTGCTGACAGGGGAGCTGAATATGTATCTGCAAGAGAGTGGATGAGGATTTGCTTTGAGCTTTTAGAGCTCTTAAAAGCCCACAAAAAGGCTATTCGTAGAGCCACAGTCAACACATTTGGTTATATTGCAAAGGCCATTGGCCCTCATGATGTATTGGCTACACTTCTGAACAACCTCAAAGTTCAAGAAAGGCAGAACAGAGTTTGTACCACTGTAGCAATAGCTATTGTTGCAGAAACATGTTCACCCTTTACAGTACTCCCTGCCTTAATGAATGAATACAGAGTTCCTGAACTGAATGTTCAAAATGGAGTGTTAAAATCGCTTTCCTTCTTGTTTGAATATATTGGTGAAATGGGAAAAGACTACATTTATGCCGTAACACCGTTACTTGAAGATGCTTTAATGGATAGAGACCTTGTACACAGACAGACGGCTAGTGCAGTGGTACAGCACATGTCACTTGGGGTTTATGGATTTGGTTGTGAAGATTCGCTGAATCACTTGTTGAACTATGTATGGCCCAATGTATTTGAGACATCTCCTCATGTAATTCAGGCAGTTATGGGAGCCCTAGAGGGCCTGAGAGTTGCTATTGGACCATGTAGAATGTTGCAATATTGTTTACAGGGTCTGTTTCACCCAGCCCGGAAAGTCAGAGATGTATATTGGAAAATTTACAACTCCATCTACATTGGTTCCCAGGACGCTCTCATAGCACATTACCCAAGAATCTACAACGATGATAAGAACACCTATATTCGTTATGAACTTGACTATATCTTATAA"
SF3B1_CDS <- s2c(SF3B1_CDS)

get_codon_from_sequence <- function(current_letter, pos, cds) {
  # This function assumes that the frame starts from position 1. 
  pos <- as.integer(pos)
  modulo_pos <- pos%%3
  
  # First position, will get the second and third letter.
  if (modulo_pos == 1){
    first_letter <- current_letter
    second_letter <- cds[pos+1]
    third_letter <- cds[pos+2]
  }
  if (modulo_pos == 2){
    first_letter <- cds[pos-1]
    second_letter <- current_letter
    third_letter <- cds[pos+1]
  }
  if (modulo_pos == 0){
    first_letter <- cds[pos-2]
    second_letter <- cds[pos-1]
    third_letter <- current_letter
  }
  
  codon <- paste0(first_letter, second_letter, third_letter)
  
  amino_acid <- GENETIC_CODE[[codon]]
  
  return(amino_acid)
}


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


### Read in dir. This is the 3rd try.
# MEK1_dir <- "/Volumes/broad_thechenlab/NextSeqOutput_Xerneas/230803_VH01286_30_AACTVW5M5/Data/Intensities/BaseCalls/"
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="CDS_pileup.tsv", full.names = T)

# This is the 2nd try.
MEK1_dir <- "/Volumes/broad_thechenlab/NextSeqOutput_Xerneas/230722_VH00997_120_AACJK37M5/Data/Intensities/BaseCalls/SF3B1/"
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="with_sequence_SF3B1_CDS_pileup.tsv", full.names = T)

# The 4th try that didn't work.
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt09_SF3B1_HEK_flow_sort_try4/"
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)
# MEK1_dir <- "/Volumes/broad_thechenlab/NextSeqOutput_Xerneas/230722_VH00997_120_AACJK37M5/Data/Intensities/BaseCalls/SF3B1/"


# We will process the 2nd try data first. 
MEK1_pileup <- read_tsv(MEK1_pileup_filenames[1]) %>% 
  mutate(sample = str_extract(basename(MEK1_pileup_filenames[1]), ".+(S\\d+)")) %>% 
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))

for (i in 2:length(MEK1_pileup_filenames)){
  tmp <- read_tsv(MEK1_pileup_filenames[i]) %>% 
    mutate(sample = str_extract(basename(MEK1_pileup_filenames[i]), ".+(S\\d+)")) %>% 
    mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
  MEK1_pileup <- rbind(MEK1_pileup, tmp)
}

# MEK1_CtoT_rep_average <- get_CtoT_average(MEK1_pileup)
# MEK1_GtoA_rep_average <- get_GtoA_average(MEK1_pileup)

MEK1_high_coverage <- MEK1_pileup %>% 
  filter((A+C+G+`T`)>10000) %>% 
  mutate(base_sum = (A+C+G+`T`)) %>%
  mutate(A_count = A) %>% 
  mutate(C_count = C) %>% 
  mutate(G_count = G) %>% 
  mutate(T_count = `T`) %>% 
  mutate(A = A/base_sum*10e6) %>% 
  mutate(C = C/base_sum*10e6) %>% 
  mutate(G = G/base_sum*10e6) %>% 
  mutate(`T` = `T`/base_sum*10e6) %>%
  mutate(base_sum = (A+C+G+`T`))


ref_all <- MEK1_high_coverage %>% 
  mutate(A_perc = A/base_sum) %>%
  mutate(C_perc = C/base_sum) %>%
  mutate(G_perc = G/base_sum) %>% 
  mutate(T_perc = `T`/base_sum) %>% 
  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
         C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
         G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
         T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum))

# Get the ref/alt percentages. 
As <- ref_all %>% filter(ref_base == "A") %>% 
  mutate(Allele_ref_count = A_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + T_count)
Cs <- ref_all %>% filter(ref_base == "C") %>% 
  mutate(Allele_ref_count = C_count) %>% 
  mutate(Allele_alt_count = A_count + G_count + T_count)
Gs <- ref_all %>% filter(ref_base == "G") %>% 
  mutate(Allele_ref_count = G_count) %>% 
  mutate(Allele_alt_count = C_count + A_count + T_count)
Ts <- ref_all %>% filter(ref_base == "T") %>% 
  mutate(Allele_ref_count = T_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + A_count)
out_dataframe_filtered <- rbind(As, Cs, Gs, Ts) %>% 
  mutate(rep = str_extract(condition, "R1|R2|R3|R4")) %>% 
  mutate(sort_group = str_extract(condition, "PcrAminus|PcrAplus")) 

ref_percentages <- out_dataframe_filtered %>% filter(grepl("PcrAminus", condition)) %>% ungroup() %>% 
  mutate(A_ref_perc = A_perc) %>% 
  mutate(C_ref_perc = C_perc) %>% 
  mutate(G_ref_perc = G_perc) %>% 
  mutate(T_ref_perc = T_perc) %>% 
  mutate(A_ref_count = A_count) %>% 
  mutate(C_ref_count = C_count) %>% 
  mutate(G_ref_count = G_count) %>% 
  mutate(T_ref_count = T_count) %>% 
  mutate(REF_Allele_ref_count = Allele_ref_count) %>% 
  mutate(REF_Allele_alt_count = Allele_alt_count) %>% 
  select(n_base, sort_group, rep, A_ref_count, C_ref_count, G_ref_count, T_ref_count, A_ref_perc, C_ref_perc, G_ref_perc, T_ref_perc, REF_Allele_ref_count, REF_Allele_alt_count) 

treated_conditions_only <- out_dataframe_filtered %>% 
  filter(grepl("PcrAplus", condition)) %>% 
  select(ref_base, condition, sort_group, rep, n_base, ref_base, A_count, C_count, G_count, T_count, A_perc, C_perc, G_perc, T_perc, Allele_ref_count, Allele_alt_count) 


treated_conditions_merged_withMinus <- merge(treated_conditions_only, ref_percentages, by = c("n_base", "rep")) %>% as_tibble() %>% arrange(rep, n_base)
# Look at mutation rate changes: 
with_ref_alt_mut_rates <- treated_conditions_merged_withMinus %>% 
  mutate(treated_count_mut_rate = Allele_alt_count/(Allele_alt_count + Allele_ref_count)) %>% 
  mutate(control_count_mut_rate = REF_Allele_alt_count/(REF_Allele_alt_count + REF_Allele_ref_count)) 

# Do fishers statistical testing.
with_ref_alt_mut_rates$fisher_res <- apply(with_ref_alt_mut_rates, 1, function(x){
  p1 <- x[['REF_Allele_ref_count']]
  p2 <- x[['Allele_ref_count']]
  p3 <- x[['REF_Allele_alt_count']]
  p4 <- x[['Allele_alt_count']]
  # print(paste(p1,p2,p3,p4))
  if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
    return(NA)
  }
  
  dat <- data.frame(
    "Control" = c(as.integer(p1), as.integer(p3)),
    "Treated" = c(as.integer(p2), as.integer(p4)),
    row.names = c("ref", "alt"),
    stringsAsFactors = FALSE
  )
  # print(dat)
  # mat<- matrix(unlist(dat), 2)
  # print(mat)
  test <- fisher.test(dat)
  return(test$p.value)
})
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% 
  mutate(fisher_res.adj = p.adjust(fisher_res, method = "BH", n = nrow(with_ref_alt_mut_rates))) %>% 
  arrange(desc(fisher_res.adj)) %>% 
  mutate(treated_to_control = treated_count_mut_rate/control_count_mut_rate)

# Add Alternate base
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% group_by(rep, n_base) %>% rowwise() %>%
  mutate(alt_base = {x <- c_across(c('A_count', 'C_count','G_count', 'T_count'));
  if (sum(!is.na(x)) >= 2) tail(head(c('A_count', 'C_count','G_count', 'T_count')[order(x, decreasing = T)],2),1) else NA})

with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% mutate(alt_base = gsub("_count", "", alt_base))
with_ref_alt_mut_rates$ref_aa <- apply(with_ref_alt_mut_rates, 1, function(x){get_codon_from_sequence(x[['ref_base']], x[['n_base']], SF3B1_CDS)})
with_ref_alt_mut_rates$alt_aa <- apply(with_ref_alt_mut_rates, 1, function(x){get_codon_from_sequence(x[['alt_base']], x[['n_base']], SF3B1_CDS)})
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% 
  mutate(mut_type = ifelse(ref_aa == alt_aa, "Silent", "Missense")) %>%
  mutate(aa_num = as.integer((n_base + 2)/3))


with_fisherres_pivoted_by_reps <- with_ref_alt_mut_rates %>% 
  select(n_base, rep, ref_base, alt_base, ref_aa, alt_aa, aa_num, mut_type, treated_count_mut_rate, control_count_mut_rate, fisher_res.adj) %>% 
  pivot_wider(names_from = c("rep"), values_from = c("treated_count_mut_rate", "control_count_mut_rate", "fisher_res.adj"))


with_fisherres_pivoted_by_reps <- with_fisherres_pivoted_by_reps %>%
  mutate(treated_mut_average = (treated_count_mut_rate_R1 + treated_count_mut_rate_R2 + treated_count_mut_rate_R3)/4) %>%
  mutate(control_mut_average = (control_count_mut_rate_R1 + control_count_mut_rate_R2 + control_count_mut_rate_R3)/4) %>% 
  mutate(treated_to_control = treated_mut_average/control_mut_average) %>% 
  mutate(log_fold_change = log(treated_to_control))

# SAVE THIS. 
try2_averaged <- with_fisherres_pivoted_by_reps %>% filter(n_base >= 1620 & n_base < 2370)


ggplot(with_fisherres_pivoted_by_reps %>% filter(n_base >= 1720 & n_base < 2370), 
       aes(n_base, log_fold_change, color = alt_base, label = aa_num)) + 
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_hline(yintercept = 2, linetype='dotted') +
  geom_text()

newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")

write.csv(with_fisherres_pivoted_by_reps, file.path(out_dir, "fig5c_SF3B1_R2_plotted_points.csv"))

g1 <- ggplot(with_fisherres_pivoted_by_reps %>% filter(n_base >= 1620 & n_base < 2370), 
       aes(n_base, log_fold_change, color = alt_base, label = aa_num)) + 
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_hline(yintercept = 2, linetype='dotted') +
  geom_point() +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_x_continuous(limits = c(1600,2400), n.breaks = 10) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("CDS base position") + 
  ylab("log(GFP+/GFP-)") 
ggsave(file.path(out_dir, "fig5c_SF3B1_wide.pdf"), g1, units = "cm", width = 16, height = 8, dpi = 300)

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370), 
       aes(control_count_mut_rate, treated_count_mut_rate, color = mut_type)) + 
  geom_point() + facet_wrap(~rep) + scale_x_log10() + scale_y_log10()

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370) , aes(log10(treated_to_control), -log10(fisher_res.adj))) + geom_point() +
  facet_wrap(~rep)

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370), aes(n_base, log10(treated_to_control), color = mut_type)) + geom_point() +
  facet_wrap(~rep)
# This is pivoted average.
ggplot(with_fisherres_pivoted_by_reps %>% filter(n_base > 1650 & n_base < 2370) , aes(n_base, log10(treated_to_control))) + geom_point() 
test <- with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370) %>% 
  filter(treated_to_control > 10)


test <- with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370)

####### LOOK AT THIRD TRY ########
# Third try
MEK1_dir <- "/Volumes/broad_thechenlab/NextSeqOutput_Xerneas/230803_VH01286_30_AACTVW5M5/Data/Intensities/BaseCalls/"
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="CDS_pileup.tsv", full.names = T)

MEK1_pileup <- read_tsv(MEK1_pileup_filenames[1]) %>% 
  mutate(sample = str_extract(basename(MEK1_pileup_filenames[1]), ".+(S\\d+)")) %>% 
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))

for (i in 2:length(MEK1_pileup_filenames)){
  tmp <- read_tsv(MEK1_pileup_filenames[i]) %>% 
    mutate(sample = str_extract(basename(MEK1_pileup_filenames[i]), ".+(S\\d+)")) %>% 
    mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)"))
  MEK1_pileup <- rbind(MEK1_pileup, tmp)
}

# MEK1_CtoT_rep_average <- get_CtoT_average(MEK1_pileup)
# MEK1_GtoA_rep_average <- get_GtoA_average(MEK1_pileup)

MEK1_high_coverage <- MEK1_pileup %>% 
  filter((A+C+G+`T`)>10000) %>% 
  mutate(base_sum = (A+C+G+`T`)) %>%
  mutate(A_count = A) %>% 
  mutate(C_count = C) %>% 
  mutate(G_count = G) %>% 
  mutate(T_count = `T`) %>% 
  mutate(A = A/base_sum*10e6) %>% 
  mutate(C = C/base_sum*10e6) %>% 
  mutate(G = G/base_sum*10e6) %>% 
  mutate(`T` = `T`/base_sum*10e6) %>%
  mutate(base_sum = (A+C+G+`T`)) %>% 
  filter(grepl("pool", sample))


ref_all <- MEK1_high_coverage %>% 
  mutate(A_perc = A/base_sum) %>%
  mutate(C_perc = C/base_sum) %>%
  mutate(G_perc = G/base_sum) %>% 
  mutate(T_perc = `T`/base_sum) %>% 
  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
         C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
         G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
         T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum))

# Get the ref/alt percentages. 
As <- ref_all %>% filter(ref_base == "A") %>% 
  mutate(Allele_ref_count = A_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + T_count)
Cs <- ref_all %>% filter(ref_base == "C") %>% 
  mutate(Allele_ref_count = C_count) %>% 
  mutate(Allele_alt_count = A_count + G_count + T_count)
Gs <- ref_all %>% filter(ref_base == "G") %>% 
  mutate(Allele_ref_count = G_count) %>% 
  mutate(Allele_alt_count = C_count + A_count + T_count)
Ts <- ref_all %>% filter(ref_base == "T") %>% 
  mutate(Allele_ref_count = T_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + A_count)
out_dataframe_filtered <- rbind(As, Cs, Gs, Ts) %>% 
  mutate(rep = str_extract(condition, "R1|R2|R3|R4")) %>% 
  mutate(sort_group = str_extract(condition, "PcrAminus|PcrAplus")) 

ref_percentages <- out_dataframe_filtered %>% filter(grepl("PcrAminus", condition)) %>% ungroup() %>% 
  mutate(A_ref_perc = A_perc) %>% 
  mutate(C_ref_perc = C_perc) %>% 
  mutate(G_ref_perc = G_perc) %>% 
  mutate(T_ref_perc = T_perc) %>% 
  mutate(A_ref_count = A_count) %>% 
  mutate(C_ref_count = C_count) %>% 
  mutate(G_ref_count = G_count) %>% 
  mutate(T_ref_count = T_count) %>% 
  mutate(REF_Allele_ref_count = Allele_ref_count) %>% 
  mutate(REF_Allele_alt_count = Allele_alt_count) %>% 
  select(n_base, sort_group, rep, A_ref_count, C_ref_count, G_ref_count, T_ref_count, A_ref_perc, C_ref_perc, G_ref_perc, T_ref_perc, REF_Allele_ref_count, REF_Allele_alt_count) 

treated_conditions_only <- out_dataframe_filtered %>% 
  filter(grepl("PcrAplus", condition)) %>% 
  select(ref_base, condition, sort_group, rep, n_base, ref_base, A_count, C_count, G_count, T_count, A_perc, C_perc, G_perc, T_perc, Allele_ref_count, Allele_alt_count) 


treated_conditions_merged_withMinus <- merge(treated_conditions_only, ref_percentages, by = c("n_base", "rep")) %>% as_tibble() %>% arrange(rep, n_base)
# Look at mutation rate changes: 
with_ref_alt_mut_rates <- treated_conditions_merged_withMinus %>% 
  mutate(treated_count_mut_rate = Allele_alt_count/(Allele_alt_count + Allele_ref_count)) %>% 
  mutate(control_count_mut_rate = REF_Allele_alt_count/(REF_Allele_alt_count + REF_Allele_ref_count)) 

# Do fishers statistical testing.
with_ref_alt_mut_rates$fisher_res <- apply(with_ref_alt_mut_rates, 1, function(x){
  p1 <- x[['REF_Allele_ref_count']]
  p2 <- x[['Allele_ref_count']]
  p3 <- x[['REF_Allele_alt_count']]
  p4 <- x[['Allele_alt_count']]
  # print(paste(p1,p2,p3,p4))
  if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
    return(NA)
  }
  
  dat <- data.frame(
    "Control" = c(as.integer(p1), as.integer(p3)),
    "Treated" = c(as.integer(p2), as.integer(p4)),
    row.names = c("ref", "alt"),
    stringsAsFactors = FALSE
  )
  # print(dat)
  # mat<- matrix(unlist(dat), 2)
  # print(mat)
  test <- fisher.test(dat)
  return(test$p.value)
})
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% 
  mutate(fisher_res.adj = p.adjust(fisher_res, method = "BH", n = nrow(with_ref_alt_mut_rates))) %>% 
  arrange(desc(fisher_res.adj)) %>% 
  mutate(treated_to_control = treated_count_mut_rate/control_count_mut_rate)

# Add Alternate base
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% group_by(rep, n_base) %>% rowwise() %>%
  mutate(alt_base = {x <- c_across(c('A_count', 'C_count','G_count', 'T_count'));
  if (sum(!is.na(x)) >= 2) tail(head(c('A_count', 'C_count','G_count', 'T_count')[order(x, decreasing = T)],2),1) else NA})

with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% mutate(alt_base = gsub("_count", "", alt_base))
with_ref_alt_mut_rates$ref_aa <- apply(with_ref_alt_mut_rates, 1, function(x){get_codon_from_sequence(x[['ref_base']], x[['n_base']], SF3B1_CDS)})
with_ref_alt_mut_rates$alt_aa <- apply(with_ref_alt_mut_rates, 1, function(x){get_codon_from_sequence(x[['alt_base']], x[['n_base']], SF3B1_CDS)})
with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% 
  mutate(mut_type = ifelse(ref_aa == alt_aa, "Silent", "Missense")) %>%
  mutate(aa_num = as.integer((n_base + 2)/3))

with_fisherres_pivoted_by_reps <- with_ref_alt_mut_rates %>% 
  select(n_base, rep, ref_base, alt_base, ref_aa, alt_aa, aa_num, mut_type, treated_count_mut_rate, control_count_mut_rate, fisher_res.adj) %>% 
  pivot_wider(names_from = c("rep"), values_from = c("treated_count_mut_rate", "control_count_mut_rate", "fisher_res.adj"))

with_fisherres_pivoted_by_reps <- with_fisherres_pivoted_by_reps %>%
  mutate(treated_mut_average = (treated_count_mut_rate_R1 + treated_count_mut_rate_R2 + treated_count_mut_rate_R3)/4) %>%
  mutate(control_mut_average = (control_count_mut_rate_R1 + control_count_mut_rate_R2 + control_count_mut_rate_R3)/4) %>% 
  mutate(treated_to_control = treated_mut_average/control_mut_average) %>% 
  mutate(log_fold_change = log(treated_to_control))

# SAVE THIS. 
try3_averaged <- with_fisherres_pivoted_by_reps %>% filter(n_base >= 1640 & n_base < 2370)
p2 <- ggplot(try2_averaged, 
             aes(n_base, log_fold_change, color = alt_base, label = aa_num)) + 
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_hline(yintercept = 2, linetype='dotted') +
  geom_text()
p3 <- ggplot(try3_averaged, 
       aes(n_base, log_fold_change, color = alt_base, label = aa_num)) + 
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_hline(yintercept = 2, linetype='dotted') +
  geom_text()


try2and3merged <- merge(try2_averaged, try3_averaged, by = "n_base")
highlight_data2 <- try2and3merged %>% filter(n_base %in% c(1849,1851,1868, 1996, 1682)) 

g2 <- ggplot(try2and3merged, aes(log_fold_change.x, log_fold_change.y, label = n_base)) + 
  geom_point(color = "gray") + geom_label() + 
  geom_point(data = highlight_data2, shape = 21, size = 1.75, color = "#e64c35", fill = "#e64c35") +
  geom_text_repel(data = highlight_data2, aes(label = n_base)) + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("log(GFP+/GFP-) Replicate 1") + 
  ylab("log(GFP+/GFP-) Replicate 2") 
ggsave(file.path(out_dir, "fig5c_SF3B1_replicate_gray.pdf"), g2, units = "cm", width = 12, height = 11, dpi = 300)



ggplot(with_fisherres_pivoted_by_reps %>% filter(n_base > 1650 & n_base < 2370), 
       aes(n_base, log10(treated_to_control), color = mut_type, label = aa_num))  + geom_text()

ggplot(with_ref_alt_mut_rates , 
       aes(control_count_mut_rate, treated_count_mut_rate, color = fisher_res.adj)) + 
  geom_point() + facet_wrap(~rep) + scale_x_log10() + scale_y_log10()

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370), 
       aes(log10(treated_to_control), -log10(fisher_res.adj), color = mut_type)) + geom_point() +
  facet_wrap(~rep)

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370) , aes(n_base, log10(treated_to_control), color = mut_type)) + geom_point() +
  facet_wrap(~rep)

ggplot(with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370) , aes(n_base, log10(treated_to_control), color = fisher_res.adj)) + geom_point() +
  facet_wrap(~rep)
# This is pivoted average.
ggplot(with_fisherres_pivoted_by_reps %>% filter(n_base > 1650 & n_base < 2370) , aes(n_base, log10(treated_to_control))) + geom_point() 
test <- with_ref_alt_mut_rates %>% filter(n_base > 1650 & n_base < 2370) %>% 
  filter(treated_to_control > 10)





##### Cursorily look at cosmic results 
SF3B1_cosmic <- read_csv("~/Downloads/Gene_mutationsMon Aug  7 19_29_36 2023.csv")
cosmic_subs <- SF3B1_cosmic  %>% arrange(desc(Count)) %>% 
  mutate(aa = str_extract(`AA Mutation`, "\\d+")) %>% 
  mutate(ref_aa = str_extract(`AA Mutation`, "[A-Z](?=\\d)"))
cosmic_subs_by_position <- cosmic_subs %>% group_by(ref_aa, aa) %>% summarise(Count = sum(Count)) %>% arrange(desc(Count)) 
# Only look at mutations with at least 3 instances?
high_cosmic_subs <- cosmic_subs_by_position %>% filter(Count>=3)

high_ratio <- with_fisherres_pivoted_by_reps %>% filter(n_base >= 1620 & n_base < 2370) %>% filter(log_fold_change > 1.5)
bases_in_cosmic <- cosmic_subs %>% filter(aa %in% high_ratio$aa_num)
with_fisherres_color_by_cosmic <- with_fisherres_pivoted_by_reps %>% filter(n_base >= 1620 & n_base < 2370) %>% 
  mutate(isInCosmic = ifelse(aa_num%in% high_cosmic_subs$aa, "Clinical", "Non-clinical")) 
with_fisherres_color_by_cosmic$isInCosmic[which(with_fisherres_color_by_cosmic$log_fold_change < 1)] <- "Threshold"
with_fisherres_color_by_cosmic <- with_fisherres_color_by_cosmic %>% mutate(isInCosmic = factor(isInCosmic, levels = c("Clinical", "Non-clinical", "Threshold")))

length(unique(bases_in_cosmic$aa))
newPalette400_cosmic <- c("#F87171","#60A5FA", "gray", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")

g3 <- ggplot(with_fisherres_color_by_cosmic, 
             aes(n_base, log_fold_change, color = isInCosmic, label = aa_num)) + 
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_hline(yintercept = 2, linetype='dotted') +
  geom_point() +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_x_continuous(limits = c(1600,2400), n.breaks = 10) + 
  scale_color_manual(values = newPalette400_cosmic) + 
  xlab("CDS base position") + 
  ylab("log(GFP+/GFP-)") 
ggsave(file.path(out_dir, "fig5c_SF3B1_wide_labelled_by_cosmic.pdf"), g3, units = "cm", width = 17, height = 8, dpi = 300)
write.csv(with_fisherres_color_by_cosmic, file.path(out_dir, "fig5c_SF3B1_R2_plotted_points.csv"))











# A_to_G_only <- MEK1_high_coverage %>% 
#   filter(ref_base == "A") %>%
#   select(n_base, A, G, condition) %>% 
#   # filter(grepl("pool", condition)) %>% 
#   # filter(grepl("RNA", condition)) %>% 
#   mutate(rep = str_extract(condition, "R1|R2|R3|R4")) %>% 
#   mutate(sort_group = str_extract(condition, "PcrAminus|PcrAplus")) %>% 
#   select(-condition)
# 
# A_to_G_pivot <- A_to_G_only %>% pivot_wider(names_from = c("sort_group", "rep"), values_from = c("G", "A"), values_fn = mean) %>% 
#   unite("G_PcrAminus", starts_with("G_PcrAminus"), sep=",") %>%
#   unite("G_PcrAplus", starts_with("G_PcrAplus"), sep=",") %>%
#   unite("A_PcrAminus", starts_with("A_PcrAminus"), sep=",") %>%
#   unite("A_PcrAplus", starts_with("A_PcrAplus"), sep=",") %>%
#   unite("G_NA", starts_with("G_NA"), sep=",") %>%
#   unite("A_NA", starts_with("A_NA"), sep=",") 
# 
# 
# 
# A_to_G_pivot_only <- A_to_G_only %>% pivot_wider(names_from = c("sort_group", "rep"), values_from = c("G", "A"), values_fn = mean)  
# A_to_G_pivot_only$fisher_res_R1 <- apply(A_to_G_pivot_only, 1, function(x){
#   p1 <- x[['G_PcrAminus_R1']]
#   p2 <- x[['G_PcrAplus_R1']]
#   p3 <- x[['A_PcrAminus_R1']]
#   p4 <- x[['A_PcrAplus_R1']]
#   if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
#     return(NA)
#   }
#   
#   dat <- data.frame(
#     "PcrAminus" = c(p1, p3),
#     "PcrAplus" = c(p2, p4),
#     row.names = c("G", "A"),
#     stringsAsFactors = FALSE
#   )
#   test <- fisher.test(dat)
#   return(test$p.value)
#   
# })
# 
# A_to_G_pivot_only$fisher_res_R2 <- apply(A_to_G_pivot_only, 1, function(x){
#   p1 <- x[['G_PcrAminus_R2']]
#   p2 <- x[['G_PcrAplus_R2']]
#   p3 <- x[['A_PcrAminus_R2']]
#   p4 <- x[['A_PcrAplus_R2']]
#   if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
#     return(NA)
#   }
#   
#   dat <- data.frame(
#     "PcrAminus" = c(p1, p3),
#     "PcrAplus" = c(p2, p4),
#     row.names = c("G", "A"),
#     stringsAsFactors = FALSE
#   )
#   test <- fisher.test(dat)
#   return(test$p.value)
#   
# })
# A_to_G_pivot_only$fisher_res_R3 <- apply(A_to_G_pivot_only, 1, function(x){
#   p1 <- x[['G_PcrAminus_R3']]
#   p2 <- x[['G_PcrAplus_R3']]
#   p3 <- x[['A_PcrAminus_R3']]
#   p4 <- x[['A_PcrAplus_R3']]
#   if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
#     return(NA)
#   }
#   
#   dat <- data.frame(
#     "PcrAminus" = c(p1, p3),
#     "PcrAplus" = c(p2, p4),
#     row.names = c("G", "A"),
#     stringsAsFactors = FALSE
#   )
#   test <- fisher.test(dat)
#   return(test$p.value)
#   
# })
# 
# A_to_G_pivot_only_p_adj <- A_to_G_pivot_only %>% 
#   mutate(fisher_res_R1.adj = p.adjust(fisher_res_R1, method = "BH", n = 3*nrow(A_to_G_pivot_only))) %>%  
#   mutate(fisher_res_R2.adj = p.adjust(fisher_res_R2, method = "BH", n = 3*nrow(A_to_G_pivot_only))) %>%  
#   mutate(fisher_res_R3.adj = p.adjust(fisher_res_R3, method = "BH", n = 3*nrow(A_to_G_pivot_only))) %>% 
#   mutate(R1_ratio = G_PcrAplus_R1/(G_PcrAplus_R1 + A_PcrAplus_R1)) %>% 
#   mutate(R2_ratio = G_PcrAplus_R2/(G_PcrAplus_R2 + A_PcrAplus_R2)) %>% 
#   mutate(R3_ratio = G_PcrAplus_R3/(G_PcrAplus_R3 + A_PcrAplus_R3)) 
# 
# A_to_G_consistent <- A_to_G_pivot_only_p_adj %>% filter(fisher_res_R1.adj < 0.01 & fisher_res_R2.adj < 0.01 & fisher_res_R3.adj < 0.01) %>% 
#   mutate(avg_ratio = (R1_ratio + R2_ratio + R3_ratio)/3)
# 
# # Get ref and alt base.
# A_to_G_consistent <- A_to_G_consistent %>% 
#   mutate(ref_base = "A") %>% 
#   mutate(alt_base = "G") %>% 
#   mutate(aa_num = as.integer((n_base + 2)/3))
# 
# A_to_G_consistent$ref_aa <- apply(A_to_G_consistent, 1, function(x){get_codon_from_sequence(x[['ref_base']], x[['n_base']], SF3B1_CDS)})
# A_to_G_consistent$alt_aa <- apply(A_to_G_consistent, 1, function(x){get_codon_from_sequence(x[['alt_base']], x[['n_base']], SF3B1_CDS)})
# A_to_G_consistent_non_silent <- A_to_G_consistent %>% filter(ref_aa != alt_aa) %>% filter(aa_num < 1000)
# write.csv(A_to_G_consistent_non_silent, "~/Dropbox (Harvard University)/03Helicase/data/AtoG_consistent_sort_try2.csv")
# 
# 
# ggplot(A_to_G_pivot_only_p_adj, aes(R2_ratio, R3_ratio, label = n_base)) + 
#   geom_point() + scale_x_log10() + scale_y_log10() + # geom_label() +
#   ggtitle("Replicate Ratio Correlation A>G") 
# 
# ggplot(A_to_G_pivot_only_p_adj , aes(R1_ratio, R3_ratio)) + 
#   geom_point() + scale_x_log10() + scale_y_log10() + 
#   ggtitle("Replicate Ratio Correlation")
# 
# A_to_G_mutate <- A_to_G_pivot %>%
#   select(n_base, G_PcrAplus, A_PcrAplus, G_PcrAminus, A_PcrAminus) %>%
#   mutate(n_base = factor(n_base)) %>%
#   mutate(I_len = 1) %>%
#   mutate(S_len = 1)
# pdat <- PDseDataSetFromMat(data.frame(A_to_G_mutate))
# pairadise_output <- pairadise(pdat, numCluster = 16)
# res <- data.frame(results(pairadise_output, p.adj = "BH", sig.level = 0.01)) 
# res$n_base <- rownames(res)
# # ggplot(data.frame(res), aes(p.value)) + geom_histogram()
# 
# # Merge with initial? 
# A_to_G_merge <- merge(A_to_G_consistent, res, by = "n_base") %>%  # A_to_G_consistent %>% # 
#   mutate(R1_ratio = R1_ratio *100) %>%
#   mutate(R2_ratio = R2_ratio *100) %>%
#   mutate(R3_ratio = R3_ratio *100) %>%
#   mutate(avg_ratio = (R1_ratio + R2_ratio + R3_ratio)/3) %>%
#   mutate(ref_base = "A") %>% 
#   mutate(alt_base = "G") %>% 
#   mutate(aa_num = as.integer((n_base + 2)/3))
# 
# # high_ratio <- A_to_G_merge %>% filter(avg_ratio > 0.5) 
# high_ratio <- A_to_G_consistent %>% filter(avg_ratio > 0.005)
# high_ratio$ref_aa <- apply(high_ratio, 1, function(x){get_codon_from_sequence(x[['ref_base']], x[['n_base']], SF3B1_CDS)})
# high_ratio$alt_aa <- apply(high_ratio, 1, function(x){get_codon_from_sequence(x[['alt_base']], x[['n_base']], SF3B1_CDS)})
# high_ratio <- high_ratio %>% filter(ref_aa != alt_aa) %>% filter(aa_num < 833)
# 
# 
# SF3B1_cosmic <- read_csv("~/Downloads/Gene_mutationsMon Aug  7 19_29_36 2023.csv")
# cosmic_subs <- SF3B1_cosmic  %>% arrange(desc(Count)) %>% 
#   mutate(aa = str_extract(`AA Mutation`, "\\d+")) %>% 
#   mutate(ref_aa = str_extract(`AA Mutation`, "[A-Z](?=\\d)"))
# cosmic_subs_by_position <- cosmic_subs %>% group_by(ref_aa, aa) %>% summarise(Count = sum(Count))
# 
# bases_in_cosmic <- cosmic_subs %>% filter(aa %in% high_ratio$aa_num)
# bases_in_cosmic_aa <- merge(bases_in_cosmic %>% select(-ref_aa), high_ratio, by.y = "aa_num", by.x="aa", all.y = T)
# # bases_in_cosmic_aa <- merge(bases_in_cosmic_aa, SF3B1_exon_RNA_liftover, by.x = "Position", by.y = "RNA_coord")
# 
# write_csv(bases_in_cosmic_aa, "~/Dropbox (Harvard University)/03Helicase/data/230808_SF3B1_filtered_by_cosmic_AtoG.csv")
# 
# ggplot(cosmic_subs_by_position, aes(as.integer(aa), Count)) + geom_bar(stat = "identity")
# 
# 
# 
# ##### Read try 2 and try 3 and compare. ######
# try2 <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt09_SF3B1_HEK_flow_sort_try4/try2and3tables/try2_consistent_hits_merged.csv")
# try3 <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt09_SF3B1_HEK_flow_sort_try4/try2and3tables/try3_consistent_hits_merged.csv")
# 
# try2 <- try2 %>% filter(aa_num < 833) %>% arrange(desc(avg_ratio))
# try3 <- try3 %>% filter(aa_num < 833) %>% arrange(desc(avg_ratio))
# 
# # try2_ratio <- try2 %>% select(n_base, avg_ratio, aa_num)
# # try3_ratio <- try3 %>% select(n_base, avg_ratio)
# 
# merged <- merge(try2, try3, by = "n_base") %>% as_tibble()
# write_csv(merged, "~/Dropbox (Harvard University)/03Helicase/data/Expt09_SF3B1_HEK_flow_sort_try4/try2and3tables/try2and3_merged.csv")
# 
# f1<- ggplot(merged, aes(avg_ratio.x, avg_ratio.y, label = aa_num.x)) + geom_point() + 
#   theme_bw()+ 
#   xlab("Biological Replicate 1") + 
#   ylab("Biological Replicate 2") + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) + scale_x_log10() + scale_y_log10()
# ggsave(file.path(out_dir, "fig4f_sf3B1_rep.pdf"), f1, units = "cm", width = 8, height = 8, dpi = 300)
# 
# 
# # Look at top for both try 2 and try 3. 
# try2_top <- try2[1:15,]$n_base
# try3_top <- try3[1:15,]$n_base
# top2and3_intersection <- intersect(try2_top, try3_top)
# 
# merged_filtered <- merged %>% filter(n_base %in% top2and3_intersection)
