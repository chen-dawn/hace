library(tidyverse)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)
library(vroom)
library(data.table)
library(rcartocolor)
##### Quantifying editing rates. #####

out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"


MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt21_Timecourse_10day/sam2tsv"

##### Uncomment to process files again. This takes a while. #####

# filename <- "/Users/dawnxi/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/sam2tsv//HEK3_WT_3_S290_aligned_to_helicase_all_amplicons_sam2tsv_filtered.tsv"

get_num_substitution_per_read <- function(filename, CHROM){
  print(filename)
  sample_name <- basename(filename)
  sample_name <- str_remove(sample_name, "_aligned_to_helicase_all_amplicons_sam2tsv_filtered.tsv")
  single_file <- fread(filename)
  
  single_file_filtered <- single_file %>% 
    filter(CHROM == CHROM) %>% 
    filter(`READ-BASE` != "N") %>% 
    filter(`REF-BASE` != "N") %>%
    filter(`READ-BASE` != `REF-BASE`) %>% 
    mutate(`REF-POS1` = as.numeric(`REF-POS1`)) 
    # filter(`REF-POS1` > 45 & `REF-POS1` < 330)
  
  # If CHROM == "03_HEK3_PCR_Product", then we filter out READ-POS0 == 76. This is a SNP.
  if (CHROM == "03_HEK3_PCR_Product"){
    single_file_filtered <- single_file_filtered %>% filter(`READ-POS0` != 76)
  }
  
  num_mut_per_read <- single_file_filtered %>% group_by(CHROM, `#Read-Name`) %>% 
    summarise(mut_per_read = n()) %>% arrange(desc(mut_per_read)) %>% 
    mutate(sample = sample_name)
  
  proportions <- num_mut_per_read %>% group_by(sample, CHROM, mut_per_read) %>%
    summarise(count = n())
  return (proportions)
}

get_num_substitution_per_read_all_chrom <- function(filename){
  print(filename)
  sample_name <- basename(filename)
  sample_name <- str_remove(sample_name, "_aligned_to_helicase_all_amplicons_sam2tsv_filtered.tsv")
  single_file <- fread(filename)
  
  single_file_filtered <- single_file %>% 
    filter(`READ-BASE` != "N") %>% 
    filter(`REF-BASE` != "N") %>%
    filter(`READ-BASE` != `REF-BASE`) %>% 
    mutate(`REF-POS1` = as.numeric(`REF-POS1`)) 
  # filter(`REF-POS1` > 45 & `REF-POS1` < 330)
  
  # If CHROM == "03_HEK3_PCR_Product", then we filter out READ-POS0 == 76. This is a SNP.
  if (CHROM == "03_HEK3_PCR_Product"){
    single_file_filtered <- single_file_filtered %>% filter(`READ-POS0` != 76)
  }
  
  num_mut_per_read <- single_file_filtered %>% group_by(CHROM, `#Read-Name`) %>% 
    summarise(mut_per_read = n()) %>% arrange(desc(mut_per_read)) %>% 
    mutate(sample = sample_name)
  
  proportions <- num_mut_per_read %>% group_by(sample, CHROM, mut_per_read) %>%
    summarise(count = n())
  return (proportions)
}


# base_summary <- single_file %>% group_by(`READ-POS0`,`READ-BASE`, `REF-BASE`) %>% summarise(n =n()) %>% arrange(desc(n))
# Get num mut per read. 
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="DC562", full.names = T)
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "03_HEK3_PCR_Product")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}
num_mut_per_read_HEK3 <- num_mut_per_read

MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="DC911", full.names = T)
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "RUNX1_300bp")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}
num_mut_per_read_RUNX1 <- num_mut_per_read

MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="DC571", full.names = T)
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "10_MAP2K1_PCR_Product")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}
num_mut_per_read_MEK1 <- num_mut_per_read


MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="DC908", full.names = T)
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "11_TNF_PCR_Product")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}
num_mut_per_read_TNF <- num_mut_per_read

MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="DC909", full.names = T)
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "12_IL6_PCR_Product")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}
num_mut_per_read_IL6 <- num_mut_per_read

num_mut_per_read_all <- rbind(num_mut_per_read_HEK3, num_mut_per_read_RUNX1, num_mut_per_read_MEK1, num_mut_per_read_TNF, num_mut_per_read_IL6)

# Get num of reads per file
# reads_per_file_filenames <- list.files(file.path(MEK1_dir), pattern="fastq_num_reads", full.names = T)
# reads_per_file <- vroom(reads_per_file_filenames) %>% 
#   mutate(sample = str_extract(filename, ".*S\\d+")) %>%
#   select(-filename)
reads_per_file_filenames <- list.files(file.path(MEK1_dir), pattern="combined_idxstats", full.names = T)
reads_per_file <- vroom(reads_per_file_filenames) %>%
  mutate(sample = str_extract(Sample, ".*S\\d+")) %>%
  mutate(CHROM = Reference) %>%
  select(-Sample, -Reference)


# Merge with stats. 
num_mut_per_read <- merge(num_mut_per_read_all, reads_per_file, by = c("sample", "CHROM")) %>%
  as_tibble()

num_mut_per_read_HEK <- num_mut_per_read %>% 
  filter(grepl("DC562", sample)) %>%
  filter(CHROM == "02_HEK3_PCR_Product")

num_mut_per_read_RUNX1 <- num_mut_per_read %>%
  filter(grepl("DC911", sample)) %>%
  filter(CHROM == "RUNX1_300bp")

num_mut_per_read_MEK1 <- num_mut_per_read %>%
  filter(grepl("DC571", sample)) %>%
  filter(CHROM == "10_MAP2K1_PCR_Product")

num_mut_per_read_TNF <- num_mut_per_read %>%
  filter(grepl("DC908", sample)) %>%
  filter(CHROM == "11_TNF_PCR_Product")

num_mut_per_read_IL6 <- num_mut_per_read %>%
  filter(grepl("DC909", sample)) %>%
  filter(CHROM == "12_IL6_PCR_Product")

num_mut_per_read <- rbind(num_mut_per_read_HEK, num_mut_per_read_RUNX1, num_mut_per_read_MEK1, num_mut_per_read_TNF, num_mut_per_read_IL6)




num_mut_per_read <- num_mut_per_read %>% 
  # extract the time from D10 for example
  mutate(condition = str_extract(sample, "DC\\d+")) %>% 
  mutate(set = str_extract(sample, "Set\\d|Lenti2")) %>% 
  mutate(day = str_extract(sample, "Ctrl|D2|D4|D6|D8|D10|D11")) %>%
  mutate(loci = str_extract(sample, "DC\\d+")) %>% 
  mutate(day = ifelse(grepl("Ctrl", sample), "D0", day)) %>%
  mutate(day = factor(day, levels = c("D0","D2", "D4", "D6", "D8", "D10", "D11"))) %>%
  filter(Mapped > 10000) %>%
  mutate(read_frac = count/Mapped) %>%
  mutate(total_mut = mut_per_read * count) %>% 
  filter(!is.na(set)) %>% 
  filter(day != "D11")


mut_per_read_less_than_5 <- num_mut_per_read %>% 
  # If read has more than 5 mutations, we merge them to 6.
  mutate(mut_per_read = ifelse(mut_per_read < 6, mut_per_read, ">5")) %>% 
  # Order the factor levels.
  mutate(mut_per_read = factor(mut_per_read, levels = c("1", "2", "3", "4", "5", ">5"))) 


# Take the average of replicates
num_mut_per_read_grouped <- num_mut_per_read %>% 
  group_by(set, condition, day, mut_per_read) %>%
  summarise(mean_mut_frac = mean(read_frac)) %>%
  mutate(weighted_mut_per_read = mean_mut_frac * mut_per_read) %>% 
  ungroup() %>% 
  # If read has more than 5 mutations, we merge them to 6.
  mutate(mut_per_read = ifelse(mut_per_read < 6, mut_per_read, ">5")) %>% 
  # Order the factor levels.
  mutate(mut_per_read = factor(mut_per_read, levels = c("1", "2", "3", "4", "5", ">5")))   

newPalette400_rotate2 <- c( "#BDB2D7", "#F9BE24", "#A3E635", "#F87171","#60A5FA", "#34D399","#E879F9")
# Plot just for DC562 (HEK3).
p1 <- ggplot(num_mut_per_read_grouped %>% filter(condition == "DC562"), aes(day, mean_mut_frac, fill = mut_per_read)) + 
  geom_bar(stat = "identity", color = "black", width = 0.75, size = 0.25) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_manual(values = newPalette400_rotate2) +
  xlab("Time (Day)") + 
  ylab("Fraction of Reads")
ggsave(file.path(out_dir, "fig1g_timecourse_read_frac_HEK3.pdf"), p1, units = "cm", width = 12, height = 8, dpi = 300)


# Now we also plot for the rest
p2 <- ggplot(num_mut_per_read_grouped %>% filter(condition != "DC562"), aes(day, mean_mut_frac, fill = mut_per_read)) + 
  geom_bar(stat = "identity", color = "black", width = 0.75, size = 0.25) + 
  theme_bw() + facet_wrap(~condition, nrow = 1, scales = "free") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_manual(values = newPalette400_rotate2) +
  ylim(0,0.08) + 
  xlab("Time (Day)") + 
  ylab("Fraction of Reads")
ggsave(file.path(out_dir, "fig1g_timecourse_read_frac_others.pdf"), p2, units = "cm", width = 24, height = 6, dpi = 300)


# Calculate the mean number of reads per day.
num_mut_per_read_mean <- num_mut_per_read %>% 
  group_by(set, sample, condition, day, mut_per_read) %>%
  summarise(mean_mut_frac = mean(read_frac)) %>%
  mutate(weighted_mut_per_read = mean_mut_frac * mut_per_read) %>% 
  ungroup() %>%
  group_by(set, sample, condition, day) %>%
  summarise(mut_frac = sum(weighted_mut_per_read)) %>%
  ungroup()  %>% 
  group_by(set, condition, day) %>% 
  summarise(mean_mut_frac = mean(mut_frac), sd_mut_frac = sd(mut_frac))

p4 <- ggplot(num_mut_per_read_mean, aes(day, mean_mut_frac, group = condition, color = condition)) + 
  geom_line( size = 1) + 
  geom_point(size = 3) + 
  geom_pointrange(aes(ymin = mean_mut_frac - sd_mut_frac, ymax = mean_mut_frac + sd_mut_frac))+ 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_manual(values = newPalette400_rotate2) +
  xlab("Time (Day)") + 
  ylab("Mean Fraction of Reads")

ggsave(file.path(out_dir, "fig1g_timecourse_mean mutation reviewer.pdf"), p4, units = "cm", width = 15, height = 10, dpi = 300)


# Sum the >5 mutations
num_mut_per_read_grouped <- num_mut_per_read %>% 
  group_by(set, condition, day, mut_per_read) %>%
  summarise(mean_mut_frac = mean(read_frac)) %>%
  # mutate(weighted_mut_per_read = mean_mut_frac * mut_per_read) %>% 
  ungroup() %>% 
  # If read has more than 5 mutations, we merge them to 6.
  mutate(mut_per_read = ifelse(mut_per_read < 6, mut_per_read, ">5")) %>% 
  group_by(set, condition, day, mut_per_read) %>%
  summarise(mean_mut_frac = sum(mean_mut_frac)) %>%
  # Order the factor levels.
  mutate(mut_per_read = factor(mut_per_read, levels = c("1", "2", "3", "4", "5", ">5")))

newPalette400_rotate2 <- c( "#F9BE24", "#A3E635", "#F87171","#60A5FA", "#34D399","#E879F9")
p2 <- ggplot(num_mut_per_read_grouped %>% filter(mut_per_read != 1), aes(day, mean_mut_frac, fill = mut_per_read)) + 
  geom_bar(stat = "identity", color = "black", width = 0.75, size = 0.25) + 
  theme_bw() + facet_wrap(~condition, nrow = 1, scales = "free") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_manual(values = newPalette400_rotate2) +
ylim(0,0.04) + 
  xlab("Time (Day)") + 
  ylab("Fraction of Reads")
ggsave(file.path(out_dir, "fig1g_timecourse_read_frac_others_summed.pdf"), p2, units = "cm", width = 30, height = 6, dpi = 300)


num_mut_per_read_grouped_is1 <- num_mut_per_read_grouped %>% filter(mut_per_read == "1")
num_mut_per_read_grouped_not1 <- num_mut_per_read_grouped %>% filter(mut_per_read != "1") %>% 
  group_by(set, condition, day) %>% 
  summarise(mean_mut_frac = sum(mean_mut_frac)) %>% 
  mutate(mut_per_read = "2-6")

# Left join the 2 dataframes
num_mut_per_read_grouped_ratio<- num_mut_per_read_grouped_is1 %>% 
  left_join(num_mut_per_read_grouped_not1, by = c("set", "condition", "day")) %>% 
  mutate(mean_mut_frac_ratio = mean_mut_frac.y/mean_mut_frac.x) 

# Rename the conditons
num_mut_per_read_grouped_ratio <- num_mut_per_read_grouped_ratio %>% 
  mutate(condition = ifelse(condition == "DC562", "HEK3", condition)) %>% 
  mutate(condition = ifelse(condition == "DC571", "MAP2K1", condition)) %>% 
  mutate(condition = ifelse(condition == "DC911", "RUNX1", condition)) %>% 
  mutate(condition = ifelse(condition == "DC908", "TNF", condition)) %>% 
  mutate(condition = ifelse(condition == "DC909", "IL6", condition))
# Plot the ratio
ggplot(num_mut_per_read_grouped_ratio, aes(day, mean_mut_frac_ratio, color = condition, group = condition)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_color_manual(values = newPalette400_rotate2) +
  xlab("Time (Day)") + 
  ylab("Ratio of reads with >1 mutations to reads with 1 mutation") +
  # Change legend title
  labs(color = "Loci")
ggsave(file.path(out_dir, "fig1g_timecourse_read_frac_ratio.pdf"), width = 6, height =4, dpi = 300)
