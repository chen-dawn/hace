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
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/RNA/"
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
# write.csv(MEK1_pileup, "~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/RNA/merged_all_A375_selection_RNA_pileups.csv")

MEK1_pileup <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt04_MEK1_A375_selection/RNA/merged_all_A375_selection_RNA_pileups.csv")

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

remove_ref_base_coverage <- function(x){
  print(x)
  if (x['ref_base'] == "A"){
    x['A'] = 0
  }
  return (x)
}


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



# Select only the S and T drug and also the conditions. Drug1 is the GFP control. 
out_dataframe_filtered <- out_dataframe %>% 
  filter(condition %in% c("A375-RNA-Ctrl", "A375-RNA-2S-day14", "A375-RNA-2S-day20", "A375-RNA-2T-day14", "A375-RNA-2T-day20")) %>%
  # filter(sample != "A375-RNA-Ctrl-B07_S11") %>% 
  mutate(cds_base = n_base - 436)
# out_dataframe_filtered$exon_num <- NA
# # merged_both_modes$exon_num[which(merged_both_modes$n_base>48090 & merged_both_modes$n_base < 48400)] <- "Exon 2"
# # merged_both_modes$exon_num[which(merged_both_modes$n_base>49810 & merged_both_modes$n_base < 50100)] <- "Exon 3"
# # merged_both_modes$exon_num[which(merged_both_modes$n_base>94640 & merged_both_modes$n_base < 94960)] <- "Exon 6"
# out_dataframe_filtered$exon_num[which(out_dataframe_filtered$n_base>=48116 & out_dataframe_filtered$n_base < 48327)] <- "Exon 2"
# out_dataframe_filtered$exon_num[which(out_dataframe_filtered$n_base>=49835 & out_dataframe_filtered$n_base < 49982)] <- "Exon 3"
# out_dataframe_filtered$exon_num[which(out_dataframe_filtered$n_base>=94844 & out_dataframe_filtered$n_base < 94969)] <- "Exon 6"
# out_dataframe_filtered <- out_dataframe_filtered %>% filter(!is.na(exon_num)) 
# # Add CDS num.
# Exon2 <- out_dataframe_filtered %>% filter(exon_num == "Exon 2") %>% mutate(cds_base = n_base + 81 - 48116)
# Exon3 <- out_dataframe_filtered %>% filter(exon_num == "Exon 3") %>% mutate(cds_base = n_base + 292 - 49835)
# Exon6 <- out_dataframe_filtered %>% filter(exon_num == "Exon 6") %>% mutate(cds_base = n_base + 569 - 94844)
# out_dataframe_filtered <- rbind(Exon2, Exon3, Exon6)

# ggplot(out_dataframe_filtered, aes(n_base, A_perc*10)) + geom_point() + facet_grid(condition~exon_num, scales = "free_x")
# ggplot(out_dataframe_filtered, aes(n_base, T_perc*10)) + geom_point() + facet_grid(condition~exon_num, scales = "free_x")

# Get the mut counts for A, C, G, T separately
As <- out_dataframe_filtered %>% filter(ref_base == "A") %>% 
  mutate(Allele_ref_count = A_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + T_count)
Cs <- out_dataframe_filtered %>% filter(ref_base == "C") %>% 
  mutate(Allele_ref_count = C_count) %>% 
  mutate(Allele_alt_count = A_count + G_count + T_count)
Gs <- out_dataframe_filtered %>% filter(ref_base == "G") %>% 
  mutate(Allele_ref_count = G_count) %>% 
  mutate(Allele_alt_count = C_count + A_count + T_count)
Ts <- out_dataframe_filtered %>% filter(ref_base == "T") %>% 
  mutate(Allele_ref_count = T_count) %>% 
  mutate(Allele_alt_count = C_count + G_count + A_count)

out_dataframe_filtered <- rbind(As, Cs, Gs, Ts)

# Reference percentages. 
ref_percentages <- out_dataframe_filtered %>% filter(condition == "A375-RNA-Ctrl") %>% ungroup() %>% 
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
  select(n_base, cds_base, A_ref_count, C_ref_count, G_ref_count, T_ref_count,
         A_ref_perc, C_ref_perc, G_ref_perc, T_ref_perc, REF_Allele_ref_count, REF_Allele_alt_count) %>% 
  group_by(n_base) %>%
  summarise(cds_base = mean(cds_base),
            A_ref_count = mean(A_ref_count),
            C_ref_count = mean(C_ref_count),
            G_ref_count = mean(G_ref_count),
            T_ref_count = mean(T_ref_count),
            A_ref_perc = mean(A_ref_perc),
            C_ref_perc = mean(C_ref_perc),
            G_ref_perc = mean(G_ref_perc),
            T_ref_perc = mean(T_ref_perc),
            REF_Allele_ref_count = mean(REF_Allele_ref_count),
            REF_Allele_alt_count = mean(REF_Allele_alt_count)) %>% ungroup()

treated_conditions_only <- out_dataframe_filtered %>% select(sample, condition, n_base, ref_base, A_count, C_count, G_count, T_count, A_perc, C_perc, G_perc, T_perc, Allele_ref_count, Allele_alt_count) %>% 
  filter(condition != "A375-RNA-Ctrl")
treated_conditions_only <- merge(treated_conditions_only, ref_percentages, by = "n_base") %>% as_tibble() %>% arrange(sample, n_base)


# Look at mutation rate changes: 
with_ref_alt_mut_rates <- treated_conditions_only %>% 
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

with_ref_alt_mut_rates <- with_ref_alt_mut_rates %>% mutate(fisher_res.adj = p.adjust(fisher_res, method = "BH", n = nrow(with_ref_alt_mut_rates)))
low_p <- with_ref_alt_mut_rates %>% filter(fisher_res.adj<0.05)

# Looks like Day 14 data is already not bad. But because this is WT control there's some skew. So we are not going to use this.
# ggplot(with_ref_alt_mut_rates, aes(control_count_mut_rate, treated_count_mut_rate, label = cds_base)) + geom_point() + # geom_label() + 
#   scale_x_log10() + scale_y_log10() + facet_wrap(~sample) 


# Only look at the day14 data and plot. 
newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")
day14s <- with_ref_alt_mut_rates %>% filter(condition =="A375-RNA-2S-day20") %>% 
  select(n_base,cds_base, A_perc, C_perc, G_perc, T_perc) %>% 
  pivot_longer(cols = contains("perc"), values_to = "mut_to_ref", names_to = "mut_mode")

g14s <- ggplot(day14s, aes(cds_base, mut_to_ref*100, fill = mut_mode)) + # geom_point() + 
  # facet_wrap(~condition, nrow = 2, scales = "free_y") + 
  # geom_point() + 
  geom_bar(stat = "identity") + # geom_hline(yintercept = 10) + 
  xlab("CDS base position") + 
  ylab("Allele Frequency (%)") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + # xlim(0, 700) + 
  scale_fill_manual(values = newPalette400_rotate) + 
  scale_x_continuous(limits = c(0,700), n.breaks = 10)


day14t <- with_ref_alt_mut_rates %>% filter(condition =="A375-RNA-2T-day20") %>% 
  select(n_base,cds_base, A_perc, C_perc, G_perc, T_perc) %>% 
  pivot_longer(cols = contains("perc"), values_to = "mut_to_ref", names_to = "mut_mode")

g14t <- ggplot(day14t, aes(cds_base, mut_to_ref*100, fill = mut_mode)) + # geom_point() + 
  # facet_wrap(~condition, nrow = 2, scales = "free_y") + 
  geom_bar(stat = "identity") + # geom_hline(yintercept = 10) + 
  xlab("CDS base position") + 
  ylab("Allele Frequency (%)") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + # xlim(0, 700) + 
  scale_fill_manual(values = newPalette400_rotate) + 
  scale_x_continuous(limits = c(0,700), n.breaks = 10)

ggsave(file.path(out_dir, "fig4a_drug_selection_allele_frequency_based_on_RNA_wide.pdf"), arrangeGrob(g14t, g14s, nrow = 2), units = "cm", width = 15, height = 10, dpi = 300)
ggsave(file.path(out_dir, "fig4a_drug_selection_allele_frequency_based_on_RNA_narrow.pdf"), arrangeGrob(g14t, g14s, nrow = 1), units = "cm", width = 22, height = 6, dpi = 300)






### Look at fold change to ref:
fold_change_to_ref <- treated_conditions_only %>%
  mutate(A_mut_to_ref = A_perc/A_ref_perc) %>% 
  mutate(C_mut_to_ref = C_perc/C_ref_perc) %>% 
  mutate(G_mut_to_ref = G_perc/G_ref_perc) %>% 
  mutate(T_mut_to_ref = T_perc/T_ref_perc) %>% 
  mutate(aa_num = as.integer((cds_base+2)/3)) %>%
  filter(!(A_mut_to_ref > 5 & A_perc < 1e-2) | !(C_mut_to_ref > 5 & C_perc < 1e-2) | (G_mut_to_ref > 5 & G_perc < 1e-2) | !(T_mut_to_ref > 5 & T_perc < 1e-2))
  
fold_change_to_ref_pivot <- fold_change_to_ref %>%
  pivot_longer(cols = contains("mut_to_ref"), values_to = "mut_to_ref", names_to = "mut_mode") 

g5_t <- ggplot(fold_change_to_ref_pivot %>% filter(condition =="mix2t-day20"), aes(cds_base, mut_to_ref, fill = mut_mode, label = ref_base)) + # geom_point() + 
  # facet_wrap(~condition, nrow = 2, scales = "free_y") + 
  geom_bar(stat = "identity") + # geom_hline(yintercept = 10) + 
  xlab("CDS base position") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 700) + 
  scale_fill_manual(values = newPalette400_rotate) + 
  scale_x_continuous(limits = c(0,700), n.breaks = 10)
g5_s <- ggplot(fold_change_to_ref_pivot %>% filter(condition =="mix2s-day20"), aes(cds_base, mut_to_ref, fill = mut_mode, label = ref_base)) + # geom_point() + 
  # facet_wrap(~condition, nrow = 2, scales = "free_y") + 
  geom_bar(stat = "identity") + # geom_hline(yintercept = 10) + 
  xlab("CDS base position") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 700) + 
  scale_fill_manual(values = newPalette400_rotate) + 
  scale_x_continuous(limits = c(0,700), n.breaks = 10)
# ggsave(file.path(out_dir, "fig4a_drug_selection.pdf"), arrangeGrob(g5_t, g5_s, nrow = 2), units = "cm", width = 11, height = 10, dpi = 300)

# filter(A_mut_to_ref > 5 & A_perc < 0.01) %>% 
  # filter(C_mut_to_ref > 5 & C_perc < 0.01) %>% 
  # filter(G_mut_to_ref > 5 & G_perc < 0.01) %>% 
  # filter(T_mut_to_ref > 5 & T_perc < 0.01) 
  
ggplot(fold_change_to_ref_pivot %>% filter(condition =="mix2s-day20"), aes(cds_base, A_perc, fill = mut_mode, label = ref_base)) + # geom_point() + 
  # facet_wrap(~condition, nrow = 2, scales = "free_y") + 
  geom_bar(stat = "identity") + # geom_hline(yintercept = 10) + 
  xlab("CDS base position") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 700) + 
  scale_fill_manual(values = newPalette400_rotate) + 
  scale_x_continuous(limits = c(0,700), n.breaks = 10)
  
ggplot(fold_change_to_ref, aes(A_perc, A_mut_to_ref, label = aa_num)) + 
  geom_point() + 
  facet_grid(condition~exon_num) + 
  geom_point(aes(C_perc, C_mut_to_ref), color = "red") + 
  geom_point(aes(G_perc, G_mut_to_ref), color = "blue") + 
  geom_point(aes(T_perc, T_mut_to_ref), color = "orange") + 
  scale_y_log10() + scale_x_log10() + ylab("Mut to Ref Fold Change") + xlab("Allele Frequency")
