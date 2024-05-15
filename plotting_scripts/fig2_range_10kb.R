library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)
library(data.table)

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
 
MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt22_range_10kb/"
# MEK1_dir <- "/Volumes/broad_thechenlab/Dawn/NextSeq/240403_VL00297_275_AAFJH5MM5/Data/Intensities/BaseCalls/testing/"
MEK1_dir <- "/Volumes/broad_thechenlab/Dawn/NextSeq/240403_VL00297_275_AAFJH5MM5/Data/Intensities/BaseCalls/Long10kb/"


MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="Long.*pileup\\.tsv$", full.names = T)

all_files <- vroom(MEK1_pileup_filenames, id = "filename")
all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>%
  mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>%
  select(-filename)

all_files_df <- all_files_df %>%
  mutate(sample = gsub("Long10kb-", "", sample)) %>%
  mutate(condition = gsub("Long10kb-", "", condition)) %>%
  mutate(condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition)) %>%
  mutate(condition = gsub("_", "-", condition)) 
  

# fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt22_range_10kb/merged_all_tagment_pileups.csv")

# c <- reconcat_all %>% group_by(sample) %>% summarise(avg = mean(base_sum))

# all_files_df <- fread("~/Dropbox (Harvard University)/03Helicase/data/Expt22_range_10kb/merged_all_tagment_pileups.csv")
MEK1_pileup <- all_files_df %>% 
  mutate(base_sum = (A + C + `T` + G )) %>%
  filter(base_sum > 8000)
# MEK1_pileup <- MEK1_pileup %>% mutate(permname = paste0(sample, chr, n_base))

## More Locations
# DC573 44231 reverse
# DC576 2972 forward
# DC577 7846 forward
DC329 <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 10641) %>% filter(chr == "TNF_10kb") %>% mutate(loci = "TNF")
# DC329 <- flip_DNA(DC329)

DC330 <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 4006) %>% filter(chr == "IL6_10kb") %>% mutate(loci = "IL6")
# DC330 <- flip_DNA(DC330)

DC573 <- MEK1_pileup %>% filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 15329) %>% filter(chr == "HEK3_10kb") %>% mutate(loci = "HEK3")

# Need to change names.
reconcat_all_10kb <- rbind(DC329, DC330, DC573) # %>% filter(n_base > -2000 & n_base < 2000)

# Also read the files from the old samples
all_files_df <- fread("~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")
MEK1_pileup <- all_files_df %>% 
  mutate(base_sum = (A + C + `T` + G )) 
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

DC329 <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 2478) %>% filter(chr == "chr6_(NC_000006)_TNF_long_extraction") %>% mutate(loci = "TNF") %>% mutate(chr = "TNF_10kb")
DC329 <- flip_DNA(DC329)

DC330 <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 1534) %>% filter(chr == "NC_000007_-_IL6_gene") %>% mutate(loci = "IL6") %>% mutate(chr = "IL6_10kb")
DC330 <- flip_DNA(DC330)

DC573 <- MEK1_pileup %>% filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 44231) %>% filter(chr == "DC573_02_HEK3") %>% mutate(loci = "HEK3") %>% mutate(chr = "HEK3_10kb")

# Need to change names.
reconcat_all_1kb <- rbind(DC329, DC330, DC573) %>% filter(n_base < 2000)


reconcat_all_1kb %>% filter(n_base > 1250 & n_base < 1265) %>% filter(chr == "TNF_10kb") %>% arrange(sample, n_base) %>% filter(condition == "DC329-DC332")
reconcat_all_10kb %>% filter(n_base > 1250 & n_base < 1265) %>% filter(chr == "TNF_10kb") %>% arrange(sample, n_base) %>% filter(condition == "DC329-DC333") 

merged_reconcat <- rbind(reconcat_all_1kb %>% mutate(batch = "1kb"), reconcat_all_10kb%>% mutate(batch = "10kb")) %>% 
  mutate(GtoA = A/(A+G)*100) %>% 
  mutate(CtoT = `T`/(C + `T`)*100) %>% 
  mutate(GtoA = ifelse(ref_base == "G", GtoA, NA)) %>% 
  mutate(CtoT = ifelse(ref_base == "C", CtoT, NA))


# Plot the coverage map for IL6 loci for the 1kb and 10kb range.
# ggplot(reconcat_all_1kb %>% filter(chr == "NC_000007_-_IL6_gene"), aes(n_base, base_sum, color = sample)) + geom_point() + facet_grid(loci ~ .)
# ggplot(reconcat_all_10kb %>% filter(chr == "IL6_10kb"), aes(n_base, base_sum, color = sample)) + geom_point() + facet_grid(loci ~ .) + xlim(0,1250)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all_10kb %>% 
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
  mutate(isControl = ifelse(grepl("DC332", condition), 0, 1))

ggplot(with_mutation_rate, aes(n_base, GtoA, color = Cas9Display)) + 
  geom_point() + 
  facet_grid(loci~helicaseDisplay)  + ylim(0,1)

##### FIRST FOR FORWARD DIRECTION ######
# grouped_by_condition <- with_mutation_rate %>% 
#   group_by(condition, n_base, chr, ref_base, loci, helicaseDisplay, Cas9Display) %>% 
#   dplyr::summarise(GtoA = mean(GtoA), CtoT = mean(CtoT), GtoA_sd = sd(GtoA), CtoT_sd = sd(CtoT)) %>% 
#   filter(n_base >= 7) %>% 
#   filter(GtoA <10 & GtoA > 0)

grouped_by_condition <- with_mutation_rate %>% 
  filter(n_base >= 7) %>% 
  filter(GtoA <10 & GtoA > 0)


##### Get Rolling Median #####
unique_condition <- unique(grouped_by_condition$sample)
all_bases <- data.frame(n_base = seq(0, 10200))
my_data <- list()
for (i in seq_along(unique_condition)){
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(sample == unique_condition[i]) 
  condition_name = tmp_mat$sample[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(sample = condition_name) %>% 
    mutate(loci = loci_name) %>% 
    mutate(helicaseDisplay = helicaseDisplay_name) %>% 
    mutate(Cas9Display = Cas9Display_name) %>% 
    arrange(n_base)
  with_rolling_mean <- tmp_mat %>% 
    mutate(rollingGtoA = zoo::rollapply(GtoA, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) %>% 
    mutate(rollingCtoT = zoo::rollapply(CtoT, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) 
  
  my_data[[i]] <- with_rolling_mean
}
all_roll_means <- do.call("rbind", my_data)

# DC573_only <- all_roll_means %>% filter(loci == "HEK3")
# DC573_only_all_base <- grouped_by_condition %>% filter(loci == "HEK3")
# ggplot(DC573_only_all_base, aes(n_base, GtoA)) + geom_point() + facet_wrap(~condition)


##### Fit Best fit line #####
data_for_fit_forward <- all_roll_means %>% filter(!is.na(helicaseDisplay)) %>% 
  mutate(logRollingGtoA = log(rollingGtoA)) %>% 
  # Reassign condition from sample name
  mutate(condition = str_extract(sample, ".*(?=-[[:upper:]]\\d+_S\\d+)")) %>% 
  # Convert all NaN to NA
  mutate(rollingGtoA = ifelse(is.nan(rollingGtoA), NA, rollingGtoA)) 

spacing100_nbase <- seq(100,12000, by = 200)
forward_points <- data_for_fit_forward %>% filter(n_base %in% spacing100_nbase)#  %>% filter(loci == "HEK3" & grepl("DC333", condition))

forward_points_nCas9 <- forward_points %>% filter(grepl("DC332", sample))
forward_points_control <- forward_points %>% filter(!grepl("DC332", sample)) %>%
  filter(!is.na(rollingGtoA))%>%
  dplyr::group_by(n_base,condition, loci, helicaseDisplay) %>% 
  dplyr::summarise(rollingGtoA = mean(rollingGtoA, na.rm = T), n=n()) %>% ungroup()

  

t_test_one_sample <- function(x){
  if (length(x) < 2){
    return(NA)
  }
  
  return(t.test(x, mu = 0)$p.value)
}
merged_forward <- merge(forward_points_nCas9, forward_points_control, by = c("n_base", "loci", "helicaseDisplay")) %>% 
  # select(sample, condition.x, n_base, loci, helicaseDisplay, rollingGtoA.x, rollingGtoA.y) %>% 
  dplyr::rename(GtoA_nCas9 = rollingGtoA.x, GtoA_control = rollingGtoA.y, treatment_condition = condition.x) %>% 
  mutate(GtoA_diff = GtoA_nCas9 - GtoA_control) %>% 
  filter(GtoA_control < 0.2) %>%
  # if it's below 0, then set it to 0.
  mutate(GtoA_diff = ifelse(GtoA_diff < 0, 0, GtoA_diff)) %>% arrange(treatment_condition, n_base) %>% 
  select(n_base, loci, sample, treatment_condition, GtoA_diff) %>% filter(!is.na(GtoA_diff)) %>%
  # Get the mean and standard deviation of the GtoA grouped by condition
  dplyr::group_by(treatment_condition, n_base, loci) %>%
  dplyr::summarise(avg_GtoA_diff = mean(GtoA_diff, na.rm = T), rollingGtoA_sd = sd(GtoA_diff, na.rm = T), n = n(), p.value = t_test_one_sample(GtoA_diff)) %>% 
  rename(condition = treatment_condition, GtoA_diff = avg_GtoA_diff) %>% ungroup()%>%
  # FDR P-value adjustment
  mutate(p.value.adj = p.adjust(p.value, method = "bonferroni"))


# 
# ggplot(forward_points, aes(n_base, rollingGtoA, color = Cas9Display)) + geom_point() + 
#   geom_line() + 
#   facet_grid(loci~helicaseDisplay) + ggtitle("101 rolling window width, 100 bp is one point") + theme_bw()

ggplot(merged_forward, aes(n_base, GtoA_diff)) + 
  geom_pointrange(aes(ymin = GtoA_diff - rollingGtoA_sd, ymax = GtoA_diff + rollingGtoA_sd)) +
  geom_line() + 
  facet_wrap(loci~condition) + ggtitle("101 rolling window width, 100,12000,by_200, GtoA_control<0.2") + theme_bw()

spacing100_nbase
merged_forward_filtered <- merged_forward %>% filter(loci == "HEK3")
# Write to output csv.
write_csv(merged_forward_filtered, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_GtoA.csv")


################################################
# This is the CtoT direction
################################################
grouped_by_condition <- with_mutation_rate %>% 
  filter(n_base >= 7) %>% 
  filter(CtoT <10 & CtoT > 0)


##### Get Rolling Median #####
unique_condition <- unique(grouped_by_condition$sample)
all_bases <- data.frame(n_base = seq(0, 10200))
my_data <- list()
for (i in seq_along(unique_condition)){
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(sample == unique_condition[i]) 
  condition_name = tmp_mat$sample[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(sample = condition_name) %>% 
    mutate(loci = loci_name) %>% 
    mutate(helicaseDisplay = helicaseDisplay_name) %>% 
    mutate(Cas9Display = Cas9Display_name) %>% 
    arrange(n_base)
  with_rolling_mean <- tmp_mat %>% 
    mutate(rollingGtoA = zoo::rollapply(GtoA, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) %>% 
    mutate(rollingCtoT = zoo::rollapply(CtoT, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) 
  
  my_data[[i]] <- with_rolling_mean
}
all_roll_means <- do.call("rbind", my_data)





# DC573_only <- all_roll_means %>% filter(loci == "HEK3")
# DC573_only_all_base <- grouped_by_condition %>% filter(loci == "HEK3")
# ggplot(DC573_only_all_base, aes(n_base, GtoA)) + geom_point() + facet_wrap(~condition)


##### Fit Best fit line #####
data_for_fit_forward <- all_roll_means %>% filter(!is.na(helicaseDisplay)) %>% 
  mutate(logRollingCtoT = log(rollingCtoT)) %>% 
  # Reassign condition from sample name
  mutate(condition = str_extract(sample, ".*(?=-[[:upper:]]\\d+_S\\d+)")) %>% 
  # Convert all NaN to NA
  mutate(rollingCtoT = ifelse(is.nan(rollingCtoT), NA, rollingCtoT)) 

spacing100_nbase <- seq(100,12000, by = 200)
forward_points <- data_for_fit_forward %>% filter(n_base %in% spacing100_nbase)#  %>% filter(loci == "HEK3" & grepl("DC333", condition))

forward_points_nCas9 <- forward_points %>% filter(grepl("DC332", sample))
forward_points_control <- forward_points %>% filter(!grepl("DC332", sample)) %>%
  filter(!is.na(rollingCtoT))%>%
  dplyr::group_by(n_base,condition, loci, helicaseDisplay) %>% 
  dplyr::summarise(rollingCtoT = mean(rollingCtoT, na.rm = T), n=n()) %>% ungroup()


merged_forward <- merge(forward_points_nCas9, forward_points_control, by = c("n_base", "loci", "helicaseDisplay")) %>% 
  # select(sample, condition.x, n_base, loci, helicaseDisplay, rollingGtoA.x, rollingGtoA.y) %>% 
  dplyr::rename(CtoT_nCas9 = rollingCtoT.x, CtoT_control = rollingCtoT.y, treatment_condition = condition.x) %>% 
  mutate(CtoT_diff = CtoT_nCas9 - CtoT_control) %>% 
  filter(CtoT_control < 0.2) %>%
  # if it's below 0, then set it to 0.
  mutate(GtoA_diff = ifelse(CtoT_diff < 0, 0, CtoT_diff)) %>% arrange(treatment_condition, n_base) %>% 
  select(n_base, loci, sample, treatment_condition, CtoT_diff) %>% filter(!is.na(CtoT_diff)) %>%
  # Get the mean and standard deviation of the GtoA grouped by condition
  dplyr::group_by(treatment_condition, n_base, loci) %>%
  dplyr::summarise(avg_CtoT_diff = mean(CtoT_diff, na.rm = T), rollingCtoT_sd = sd(CtoT_diff, na.rm = T), n = n(), p.value = t.test(CtoT_diff, mu = 0)$p.value) %>% 
  rename(condition = treatment_condition, CtoT_diff = avg_CtoT_diff) %>% ungroup() 


# 
# ggplot(forward_points, aes(n_base, rollingGtoA, color = Cas9Display)) + geom_point() + 
#   geom_line() + 
#   facet_grid(loci~helicaseDisplay) + ggtitle("101 rolling window width, 100 bp is one point") + theme_bw()

ggplot(merged_forward, aes(n_base, CtoT_diff)) + 
  geom_pointrange(aes(ymin = CtoT_diff - rollingCtoT_sd, ymax = CtoT_diff + rollingCtoT_sd)) +
  geom_line() + 
  facet_wrap(loci~condition) + ggtitle("101 rolling window width, 100,12000,by_200, GtoA_control<0.2") + theme_bw()


merged_forward_filtered <- merged_forward %>% filter(loci == "HEK3")
# Write to output csv.
write_csv(merged_forward_filtered, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_CtoT.csv")


############ Also do this for the old data the reverse only ############

# Plot the coverage map for IL6 loci for the 1kb and 10kb range.
# ggplot(reconcat_all_1kb %>% filter(chr == "NC_000007_-_IL6_gene"), aes(n_base, base_sum, color = sample)) + geom_point() + facet_grid(loci ~ .)
# ggplot(reconcat_all_10kb %>% filter(chr == "IL6_10kb"), aes(n_base, base_sum, color = sample)) + geom_point() + facet_grid(loci ~ .) + xlim(0,1250)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all_1kb %>% 
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
  mutate(isControl = ifelse(grepl("DC332", condition), 0, 1))

ggplot(with_mutation_rate, aes(n_base, GtoA, color = Cas9Display)) + 
  geom_point() + 
  facet_grid(loci~helicaseDisplay)  + ylim(0,1)

##### FIRST FOR FORWARD DIRECTION ######
# grouped_by_condition <- with_mutation_rate %>% 
#   group_by(condition, n_base, chr, ref_base, loci, helicaseDisplay, Cas9Display) %>% 
#   dplyr::summarise(GtoA = mean(GtoA), CtoT = mean(CtoT), GtoA_sd = sd(GtoA), CtoT_sd = sd(CtoT)) %>% 
#   filter(n_base >= 7) %>% 
#   filter(GtoA <10 & GtoA > 0)

grouped_by_condition <- with_mutation_rate %>% 
  filter(n_base <= 0) %>% 
  filter(CtoT <10 & CtoT > 0)


##### Get Rolling Median #####
unique_condition <- unique(grouped_by_condition$sample)
all_bases <- data.frame(n_base = seq(-2000, 0))
my_data <- list()
for (i in seq_along(unique_condition)){
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(sample == unique_condition[i]) 
  condition_name = tmp_mat$sample[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(sample = condition_name) %>% 
    mutate(loci = loci_name) %>% 
    mutate(helicaseDisplay = helicaseDisplay_name) %>% 
    mutate(Cas9Display = Cas9Display_name) %>% 
    arrange(n_base)
  with_rolling_mean <- tmp_mat %>% 
    mutate(rollingGtoA = zoo::rollapply(GtoA, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) %>% 
    mutate(rollingCtoT = zoo::rollapply(CtoT, width = 101, mean, align = "center", na.rm = TRUE, fill = NA)) 
  
  my_data[[i]] <- with_rolling_mean
}
all_roll_means <- do.call("rbind", my_data)





# DC573_only <- all_roll_means %>% filter(loci == "HEK3")
# DC573_only_all_base <- grouped_by_condition %>% filter(loci == "HEK3")
# ggplot(DC573_only_all_base, aes(n_base, GtoA)) + geom_point() + facet_wrap(~condition)


##### Fit Best fit line #####
data_for_fit_forward <- all_roll_means %>% filter(!is.na(helicaseDisplay)) %>% 
  # Reassign condition from sample name
  mutate(condition = str_extract(sample, ".*(?=-[[:upper:]]\\d+_S\\d+)")) %>% 
  mutate(condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition)) %>% 
  # Convert all NaN to NA
  mutate(rollingCtoT = ifelse(is.nan(rollingCtoT), NA, rollingCtoT)) 

spacing100_nbase <- seq(-2000,0, by = 100)
forward_points <- data_for_fit_forward %>% filter(n_base %in% spacing100_nbase) %>% filter(loci == "HEK3")

forward_points_nCas9 <- forward_points %>% filter(grepl("DC332", sample))
forward_points_control <- forward_points %>% filter(condition %in% c("DC573_DC333", "DC573_DC334","DC573_DC335", "DC573_pBO101")) %>%
  filter(!is.na(rollingCtoT))%>%
  dplyr::group_by(n_base,condition, loci, helicaseDisplay) %>% 
  dplyr::summarise(rollingCtoT = mean(rollingCtoT, na.rm = T), n=n()) %>% ungroup()


merged_forward <- merge(forward_points_nCas9, forward_points_control, by = c("n_base", "loci", "helicaseDisplay")) %>% 
  # select(sample, condition.x, n_base, loci, helicaseDisplay, rollingGtoA.x, rollingGtoA.y) %>% 
  dplyr::rename(CtoT_nCas9 = rollingCtoT.x, CtoT_control = rollingCtoT.y, treatment_condition = condition.x) %>% 
  mutate(CtoT_diff = CtoT_nCas9 - CtoT_control) %>% 
  filter(CtoT_control < 0.2) %>%
  # if it's below 0, then set it to 0.
  mutate(CtoT_diff = ifelse(CtoT_diff < 0, 0, CtoT_diff)) %>% arrange(treatment_condition, n_base) %>% 
  select(n_base, loci, sample, treatment_condition, CtoT_diff) %>% filter(!is.na(CtoT_diff))

  # Get the mean and standard deviation of the GtoA grouped by condition
merged_forward_grouped <- merged_forward %>% dplyr::group_by(treatment_condition, n_base, loci) %>%
  dplyr::summarise(avg_CtoT_diff = mean(CtoT_diff, na.rm = T), rollingCtoT_sd = sd(CtoT_diff, na.rm = T), n = n()) %>% 
  rename(condition = treatment_condition, CtoT_diff = avg_CtoT_diff) %>% ungroup() 


# 
# ggplot(forward_points, aes(n_base, rollingGtoA, color = Cas9Display)) + geom_point() + 
#   geom_line() + 
#   facet_grid(loci~helicaseDisplay) + ggtitle("101 rolling window width, 100 bp is one point") + theme_bw()

ggplot(merged_forward_grouped, aes(n_base, CtoT_diff)) + 
  geom_pointrange(aes(ymin = CtoT_diff - rollingCtoT_sd, ymax = CtoT_diff + rollingCtoT_sd)) +
  geom_line() + 
  facet_wrap(loci~condition) + ggtitle("101 rolling window width, 100,12000,by_200, CtoT_control< 0.2") + theme_bw()


# Write to output csv.
write_csv(merged_forward_grouped, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_1kb_reverse_CtoT.csv")




################################################
grouped_by_condition <- with_mutation_rate %>% 
  filter(n_base <= 0) %>% 
  filter(GtoA < 10 & GtoA > 0)

##### Get Rolling Median #####
unique_condition <- unique(grouped_by_condition$sample)
all_bases <- data.frame(n_base = seq(-2000, 0))
my_data <- list()

for (i in seq_along(unique_condition)) {
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(sample == unique_condition[i]) 
  condition_name = tmp_mat$sample[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = TRUE) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(sample = condition_name,
           loci = loci_name,
           helicaseDisplay = helicaseDisplay_name,
           Cas9Display = Cas9Display_name) %>% 
    arrange(n_base)
  
  with_rolling_mean <- tmp_mat %>% 
    mutate(rollingGtoA = zoo::rollapply(GtoA, width = 101, mean, align = "center", na.rm = TRUE, fill = NA))
  
  my_data[[i]] <- with_rolling_mean
}

all_roll_means <- do.call("rbind", my_data)

##### Fit Best Fit Line #####
data_for_fit_forward <- all_roll_means %>% filter(!is.na(helicaseDisplay)) %>% 
  mutate(condition = str_extract(sample, ".*(?=-[[:upper:]]\\d+_S\\d+)"),
         condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition),
         rollingGtoA = ifelse(is.nan(rollingGtoA), NA, rollingGtoA))

spacing100_nbase <- seq(-2000, 0, by = 100)
forward_points <- data_for_fit_forward %>% 
  filter(n_base %in% spacing100_nbase) %>%
  filter(loci == "HEK3")

forward_points_nCas9 <- forward_points %>% filter(grepl("DC332", sample))
forward_points_control <- forward_points %>% 
  filter(condition %in% c("DC573_DC333", "DC573_DC334","DC573_DC335", "DC573_pBO101")) %>%
  filter(!is.na(rollingGtoA)) %>%
  group_by(n_base, condition, loci, helicaseDisplay) %>% 
  summarise(rollingGtoA = mean(rollingGtoA, na.rm = TRUE), n=n()) %>% ungroup()

merged_forward <- merge(forward_points_nCas9, forward_points_control, by = c("n_base", "loci", "helicaseDisplay")) %>% 
  rename(GtoA_nCas9 = rollingGtoA.x, GtoA_control = rollingGtoA.y, treatment_condition = condition.x) %>% 
  mutate(GtoA_diff = GtoA_nCas9 - GtoA_control) %>% 
  filter(GtoA_control < 0.2) %>%
  mutate(GtoA_diff = ifelse(GtoA_diff < 0, 0, GtoA_diff)) %>% 
  arrange(treatment_condition, n_base) %>% 
  select(n_base, loci, sample, treatment_condition, GtoA_diff) %>% 
  filter(!is.na(GtoA_diff))

merged_forward_grouped <- merged_forward %>% 
  group_by(treatment_condition, n_base, loci) %>%
  # Summarise the data to get mean and standard deviation
  summarise(avg_GtoA_diff = mean(GtoA_diff, na.rm = TRUE), 
            rollingGtoA_sd = sd(GtoA_diff, na.rm = TRUE), 
            n = n()) %>% 
  rename(condition = treatment_condition, GtoA_diff = avg_GtoA_diff) %>% 
  ungroup()

# Visualize the data
ggplot(merged_forward_grouped, aes(n_base, GtoA_diff)) + 
  geom_pointrange(aes(ymin = GtoA_diff - rollingGtoA_sd, ymax = GtoA_diff + rollingGtoA_sd)) +
  geom_line() + 
  facet_wrap(~loci + condition) + 
  ggtitle("101 rolling window width, 100 to 2000 by 200, GtoA control < 0.2") + 
  theme_bw()

# Save the data to a CSV file
write_csv(merged_forward_grouped, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_1kb_reverse_GtoA.csv")

###### Make some edits

library(tidyverse)

dfGtoA <- read_csv("~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_GtoA.csv")
dfCtoT <- read_csv("~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_CtoT.csv")
spacing100_nbase <- data.frame(n_base = seq(100,12000, by = 200))

unique_conditions <- unique(dfGtoA$condition)
dfGtoA_all <- data.frame()
dfCtoT_all <- data.frame()
for (c in unique_conditions){
  filtered_dfGtoA <- merge(dfGtoA %>% filter(condition == c), spacing100_nbase, all.y=T) %>% mutate(condition = c) %>% mutate(loci = "HEK3")
  filtered_dfCtoT <- merge(dfCtoT %>% filter(condition == c), spacing100_nbase, all.y=T) %>% mutate(condition = c) %>% mutate(loci = "HEK3")
  
  dfGtoA_all = rbind(dfGtoA_all, filtered_dfGtoA)
  dfCtoT_all = rbind(dfCtoT_all, filtered_dfCtoT)
}

dfGtoA <- dfGtoA_all %>% arrange(condition, n_base)
dfCtoT <- dfCtoT_all %>% arrange(condition, n_base)

write_csv(dfGtoA, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_GtoA_updated.csv")
write_csv(dfCtoT, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_10kb_merged_forward_filtered_CtoT_updated.csv")


# Want to calculate one sample t test for each row in the data frame.
one_sample_t_test <- function(mu, sd, n){
  if (is.na(mu) | is.na(sd) | is.na(n)){
    return(NA)
  }
  mu <- as.numeric(mu)
  sd <- as.numeric(sd)
  n <- as.integer(n)
  return((mu - 0)/(sd/sqrt(n)))
}
dfGtoA$tres <- apply(dfGtoA, 1, function(x) one_sample_t_test(x["GtoA_diff"], x["rollingGtoA_sd"], x["n"]))
dfCtoT$tres <- apply(dfCtoT, 1, function(x) one_sample_t_test(x["CtoT_diff"], x["rollingCtoT_sd"], x["n"]))
