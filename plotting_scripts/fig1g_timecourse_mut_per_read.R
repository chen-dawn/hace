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

MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/sam2tsv/"

##### Uncomment to process files again. This takes a while. #####
MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="HEK3", full.names = T)

filename <- "/Users/dawnxi/Dropbox (Harvard University)/03Helicase/data/Expt12_Timecourse/sam2tsv//HEK3_WT_3_S290_aligned_to_helicase_all_amplicons_sam2tsv_filtered.tsv"

get_num_substitution_per_read <- function(filename, CHROM){
  print(filename)
  sample_name <- basename(filename)
  sample_name <- str_extract(sample_name, ".*S\\d+")
  single_file <- fread(filename)
  
  single_file_filtered <- single_file %>% 
    filter(CHROM == CHROM) %>% 
    filter(`READ-BASE` != "N") %>% 
    filter(`READ-BASE` != `REF-BASE`) %>% 
    filter(`REF-POS1` > 45 & `REF-POS1` < 330)
  
  num_mut_per_read <- single_file_filtered %>% group_by(`#Read-Name`) %>% 
    summarise(mut_per_read = n()) %>% arrange(desc(mut_per_read)) %>% 
    mutate(sample = sample_name)
  
  proportions <- num_mut_per_read %>% group_by(sample, mut_per_read) %>%
    summarise(count = n())
  return (proportions)
}

# Get num mut per read. 
num_mut_per_read <- data.frame()
for (filename in MEK1_pileup_filenames){
  tmp <- get_num_substitution_per_read(filename, CHROM = "03_HEK3_PCR_Product")
  num_mut_per_read <- rbind(num_mut_per_read, tmp)
}

# Get num of reads per file
reads_per_file_filenames <- list.files(file.path(MEK1_dir), pattern="fastq_num_reads", full.names = T)
reads_per_file <- vroom(reads_per_file_filenames) %>% 
  mutate(sample = str_extract(filename, ".*S\\d+")) %>%
  select(-filename)

# Merge with stats. 
num_mut_per_read <- merge(num_mut_per_read, reads_per_file, by = "sample") %>%
  as_tibble()

num_mut_per_read <- num_mut_per_read %>% 
  mutate(rep = str_extract(sample, "\\d(?=_S\\d)")) %>% 
  mutate(time = str_extract(sample, "\\d+(?=H_)")) %>% 
  mutate(condition = str_extract(sample, ".*(?=_\\d_S\\d)")) %>% 
  mutate(time = as.integer(time)) %>%
  mutate(time = ifelse(is.na(time), 0, time)) %>% 
  filter(time < 100) %>% 
  filter(num_reads > 10000) %>%
  mutate(read_frac = count/num_reads)%>%
  mutate(total_mut = mut_per_read * count)

num_mut_per_read_grouped <- num_mut_per_read %>% 
  group_by(condition, rep, time, num_reads) %>% 
  summarise(total_mut_sum = sum(total_mut)) %>% 
  mutate(mut_per_read = total_mut_sum/num_reads)
  

num_mut_per_read_grouped_by_time <- num_mut_per_read_grouped %>% 
  group_by(condition, time) %>%
  summarise(mean_mut_per_read = mean(mut_per_read), sd_mut_per_read = sd(mut_per_read, na.rm = T)) %>%
  filter(grepl("GSPcrA", condition)| grepl("WT", condition)) %>%
  mutate(time_string = as.character(time)) 

f1 <- ggplot(num_mut_per_read_grouped_by_time, aes(time_string, mean_mut_per_read, group = 1)) + 
  geom_line(color ="#2185C5") + 
  geom_errorbar(aes(ymin=mean_mut_per_read-sd_mut_per_read, ymax=mean_mut_per_read+sd_mut_per_read), color = "#2185C5", width = 0.1) + 
  geom_point(color ="#2185C5") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("Time (h)") + 
  ylab("Number of mutations per contiguous read") 
ggsave(file.path(out_dir, "fig1g_timecourse2.pdf"), f1, units = "cm", width = 11, height = 8, dpi = 300)


### OLD DATA
# names(num_mut_per_read) <- c("ReadName", "n", "sample")
ggplot(num_mut_per_read %>% filter(mut_per_read<10) %>% 
  filter(grepl("GSPcrA", sample) | grepl("WT", sample)), aes(mut_per_read, read_frac)) + 
  geom_bar(stat = "identity") + 
  facet_grid(time~rep)

num_mut_per_read_rep_avg <- num_mut_per_read %>% 
  group_by(condition, time, mut_per_read) %>%
  summarise(read_frac = mean(read_frac)) %>% 
  mutate(mut_per_read = ifelse(mut_per_read < 6, mut_per_read, 6)) %>% 
  group_by(condition, time, mut_per_read) %>%
  summarise(read_frac = sum(read_frac)) %>%
  filter(mut_per_read > 1) %>% 
  filter(grepl("GSPcrA", condition)| grepl("WT", condition)) %>%
  mutate(time_string = as.character(time)) %>%
  mutate(mut_per_read = as.character(mut_per_read))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

newPalette400 <- c("#F87171", "#22D3EE", "#FBBF24", "#60A5FA", "#A3E635", "#A78BFA", "#34D399", "#E879F9")
newPalette500 <- c("#EF4444", "#06B6D4", "#F59E0B", "#3B82F6", "#84CC16", "#8B5CF6", "#10B981", "#D946EF")
newPalette600 <- c("#DC2626", "#0891B2", "#D97706", "#2563EB", "#65A30D", "#9333EA", "#059669", "#C026D3")
newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399", "#A78BFA", "#FBBF24", "#A3E635", "#E879F9")

library(khroma)

f1 <- ggplot(num_mut_per_read_rep_avg, aes(time_string, read_frac*100, color = mut_per_read, group = mut_per_read)) + 
  geom_point() + geom_line(size = 1.1) + 
  # geom_bar(stat = "identity", position = position_stack(reverse = TRUE), color = "black") +
  scale_color_manual(values = newPalette400_rotate) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("Time (h)") + 
  ylab("Read Percentage")

ggsave(file.path(out_dir, "fig1g_timecourse.pdf"), f1, units = "cm", width = 11, height = 8, dpi = 300)

  # theme(legend.position = c(0.9, 0.85)) 
  # khroma::scale_fill_light()
# scale_fill_brewer(palette = "YlGnBu", direction = -1)
  # scale_fill_brewer(name = "Num Mut", palette = "Set2", direction = 1)
