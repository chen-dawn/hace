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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt18_A375_BE_selection/"

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
  filter(base_sum > 1000) %>% 
  mutate(day = str_extract(condition, "D0|D14")) %>% 
  mutate(drug = str_extract(condition, "(?<=D14)\\S")) %>%
  mutate(plasmid = str_extract(condition, "DC\\d+")) 

legends <- read_csv("~/Dropbox (Harvard University)/03Helicase/data/Expt18_A375_BE_selection/A375_BE_Plasmids.csv")
MEK1_pileup <- merge(MEK1_pileup, legends, by = "plasmid")
  
with_mut_rate <- MEK1_pileup %>% mutate(base_sum = (A+C+G+`T`)) %>% 
  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
         C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
         G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
         T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum)) %>%
  mutate(mut_perc= A_perc + C_perc + G_perc + T_perc)

# ggplot(with_mut_rate, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") + scale_y_log10()

control_rates <- with_mut_rate %>% filter(day == "D0") %>% 
  ungroup() %>% 
  mutate(A_ref_perc = A_perc, C_ref_perc = C_perc, G_ref_perc = G_perc, T_ref_perc = T_perc) %>% 
  select(plasmid, n_base, A_ref_perc, C_ref_perc, G_ref_perc, T_ref_perc)  %>% 
  group_by(plasmid, n_base) %>% 
  summarise(A_ref_perc = mean(A_ref_perc), C_ref_perc = mean(C_ref_perc), G_ref_perc = mean(G_ref_perc), T_ref_perc = mean(T_ref_perc))
    

not_control <- with_mut_rate %>% filter(day == "D14")

merged_with_control <- merge(not_control, control_rates, by=c("plasmid", "n_base")) %>% as_tibble()


high_mut <- merged_with_control %>% filter(mut_perc > 0.01)
ggplot(high_mut, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") + scale_y_log10()


# Look at the locations individually. We want to see:
# DC686 (49926)
# DC687 (49926)
# DC688 (94880.2)
# DC689 (94880.2)
# DC690 (94880)
# DC691 (94882)
# DC659 (AAVS)

# OFF TARGETS
DC686 <- with_mut_rate %>% filter(plasmid == "DC686") %>% 
  filter(n_base >49920 & n_base < 49935)
ggplot(DC686, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

# GOOD 
DC687 <- with_mut_rate %>% filter(plasmid == "DC687") %>% 
  filter(n_base >49924 & n_base < 49930)
ggplot(DC687, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

DC688 <- with_mut_rate %>% filter(plasmid == "DC688") %>% 
  filter(n_base >94875 & n_base < 94885)
ggplot(DC688, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

DC689 <- with_mut_rate %>% filter(plasmid == "DC689") %>% 
  filter(n_base >94875 & n_base < 94885)
ggplot(DC689, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

# GOOD but missing a rep for T...?
DC690 <- with_mut_rate %>% filter(plasmid == "DC690") %>% 
  filter(n_base >94875 & n_base < 94885)
ggplot(DC690, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

# GOOD but also have off targets at other locations (94881, 94884)
DC691 <- with_mut_rate %>% filter(plasmid == "DC691") %>% 
  filter(n_base >94875 & n_base < 94885)
ggplot(DC691, aes(n_base, mut_perc, color = drug)) + geom_point() + facet_wrap(~plasmidName, scales = "free_x") 

merged_good_ones <- rbind(DC687, DC690, DC691)
write_csv(merged_good_ones, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_fig3j_A375_BE_validation.csv")

