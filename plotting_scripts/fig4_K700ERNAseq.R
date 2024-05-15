library(tidyverse)
library(gridExtra)
library(grid)
library(gt)
library(matrixTests)
library(ggrepel)
library(gridExtra)
library("viridis")
library(ggsci)
library(ggpubr)
library(cowplot)
out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"

rmats_df <- read_tsv("~/Dropbox (Harvard University)/03Helicase/data/Expt19_SF3B1_RMATS/A3SS.MATS.JC.txt")
# rmats_df <- read_tsv("~/Dropbox (Harvard University)/03Helicase/data/Expt19_SF3B1_RMATS/SE.MATS.JC.txt")

names(rmats_df)[1] <- "ID"

rmats_df_processed <- rmats_df %>% arrange(FDR) %>% mutate(logFDR = -log10(FDR)) %>% 
  separate(IncLevel1, into = c("IncLevel1a", "IncLevel1b", "IncLevel1c"),sep = ",") %>% 
  separate(IncLevel2, into = c("IncLevel2a", "IncLevel2b", "IncLevel2c"),sep = ",") %>%
  mutate(IncLevel1a = as.numeric(IncLevel1a)) %>% 
  mutate(IncLevel1b = as.numeric(IncLevel1b)) %>% 
  mutate(IncLevel1c = as.numeric(IncLevel1c)) %>% 
  mutate(IncLevel2a = as.numeric(IncLevel2a)) %>% 
  mutate(IncLevel2b = as.numeric(IncLevel2b)) %>% 
  mutate(IncLevel2c = as.numeric(IncLevel2c)) %>% 
  filter(!is.na(IncLevel1a) & !is.na(IncLevel1b), !is.na(IncLevel1c), !is.na(IncLevel2a), !is.na(IncLevel2b), !is.na(IncLevel2c))


rmats_df_processed <- rmats_df_processed %>% 
  mutate(IncLevel1 = (IncLevel1a + IncLevel1b + IncLevel1c)/3) %>%
  mutate(IncLevel2 = (IncLevel2a + IncLevel2b + IncLevel2c)/3) %>% 
  mutate(Level1to2 = -log10(IncLevel1/IncLevel2))


ggplot(rmats_df_processed, aes(Level1to2, logFDR)) + geom_point()

highlight_data1 <- rmats_df_processed %>% filter(FDR < 0.05 & IncLevelDifference >0.1) 
highlight_data2 <- rmats_df_processed %>% filter(FDR < 0.05 & IncLevelDifference < -0.1) 
g <- ggplot(rmats_df_processed, aes(IncLevelDifference, logFDR)) + geom_point(color = "gray") +  
  geom_point(data = highlight_data1, shape = 21, size = 1.75, color = "#e64c35", fill = "#e64c35") +
  geom_point(data = highlight_data2, shape = 21, size = 1.75, color = "#e64c35", fill = "#e64c35") +
  # geom_text_repel(data = highlight_data2, aes(label = geneSymbol)) + 
  theme_bw()+ scale_fill_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("Delta PSI") + 
  ylab("-log(FDR)")
ggsave(file.path(out_dir, "supp_5a_K700E_rmats.pdf"), g, units = "cm", width = 16, height = 12, dpi = 150)


# bc_table_all <- read_csv("~/K562_K700E/20230130_twist_library_v3_ID_barcode_ROUT.csv")

# files_dir <- "~/K562_K700E/"
# filenames <- list.files(file.path(files_dir), pattern="misspliced_analysis.csv", full.names = T)
# 
# all_insertions <- read_csv(filenames[1]) %>% 
#   mutate(sample = str_extract(basename(filenames[1]), ".+(S\\d+)")) %>% 
#   mutate(condition = str_extract(sample, "(?<=K562_)[^_]+_[^_]+")) 
# 
# for (i in 2:length(filenames)){
#   tmp <- read_csv(filenames[i]) %>% 
#     mutate(sample = str_extract(basename(filenames[i]), ".+(S\\d+)")) %>% 
#     mutate(condition = str_extract(sample, "(?<=K562_)[^_]+_[^_]+"))
#   all_insertions <- rbind(all_insertions, tmp)
# }
# 
# all_insertions <- all_insertions %>% mutate(celltype = str_extract(condition, "^([A-Z]+(?:\\d+[A-Z])?)"))
# names(all_insertions)[1] <- "cb_umi"
# all_by_umi <- all_insertions %>% separate(cb_umi, sep = "_", into = c("cb", "umi"))

# Plot upstream insertions
all_by_umi_not_zero <- all_by_umi %>% filter(upstream_insertion_length > 0 & upstream_insertion_length < 200)
ggplot(all_by_umi_not_zero %>% mutate(celltype = factor(celltype, levels = c("WT", "K700E"))), aes(x = upstream_insertion_length, y = ..density..)) + 
  # geom_histogram(fill = "4dbcd5", binwidth = 3, color = "white") + 
  geom_histogram(fill = "4dbcd5", color = "white") + 
  facet_wrap(~celltype, ncol=1) + 
  theme_classic() + xlab("Retained Intron Length") + ylab("Normalized Density") +
  theme_bw()+ scale_fill_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank()) 

# ggsave(file="~/ManuscriptFigures/21_RetainedIntron_K700E.pdf", units = "cm", width = 13, height = 12, dpi = 300) #saves g
# ggsave(file="~/ManuscriptFigures/21_RetainedIntron_K700E.png", units = "cm", width = 13, height = 12, dpi = 300) #saves g

# We will do some filtering, everything that's less than 3 upstream insertion will be considered as 0, 
# and then everything that's more than 50 will be excluded. 
all_by_umi_filtered <- all_by_umi %>% 
  select(-downstream_insertion_length, -downstream_insertion_count) %>%
  filter(upstream_insertion_length < 50) %>% 
  mutate(upstream_insertion_length = ifelse(upstream_insertion_length<3, 0, upstream_insertion_length))

# Get the proportion of different insertion lengths.
upstream_by_cb_prop <- all_by_umi_filtered %>% group_by(condition, celltype, cb, upstream_insertion_length) %>%
  summarise(count = n()) 
upstream_by_cb_prop <- upstream_by_cb_prop %>% group_by(condition, celltype, cb) %>%
  mutate(total_count = sum(count), prop = count / sum(count)) %>% ungroup() %>% arrange(prop)

# Get proportion by cell types. 
upstream_by_cb_prop_celltype <- all_by_umi_filtered %>% group_by(celltype, cb, upstream_insertion_length) %>%
  summarise(count = n()) 
upstream_by_cb_prop_celltype <- upstream_by_cb_prop_celltype %>% group_by(celltype, cb) %>%
  mutate(total_count = sum(count), prop = count / sum(count)) %>% ungroup() %>% arrange(prop)


##### RUN FISHER'S EXACT TEST. 
# Group into with and without insertion based on cell type. 
with_insertion <- all_by_umi_filtered %>% filter(upstream_insertion_length > 0) %>% 
  group_by(celltype, cb) %>% summarise(count = n()) %>% mutate(insertion_type = "with_insertion")
without_insertion <- all_by_umi_filtered %>% filter(upstream_insertion_length == 0) %>% 
  group_by(celltype, cb) %>% summarise(count = n()) %>% mutate(insertion_type = "without_insertion")

with_and_without_insertion_merged <- rbind(with_insertion, without_insertion) %>% 
  pivot_wider(names_from = c("celltype", "insertion_type"), values_from = "count", values_fill = 0)

# Perform Fisher's exact test for each row
results <- apply(with_and_without_insertion_merged[,2:5], 1, function(row) {
  contingency_table <- matrix(row, nrow = 2)
  fisher.test(contingency_table)
})

# Extract p-values from the results
with_and_without_insertion_merged$p_val <- sapply(results, function(res) res$p.value)
with_and_without_insertion_merged$log_p_val <- -log(with_and_without_insertion_merged$p_val)

# Get only the ones with upstream insertion length > 0.
upstream_by_cb_prop_celltype_with_insertion <- upstream_by_cb_prop_celltype %>%
  filter(upstream_insertion_length > 0)

# Get normalized insertion length for both WT and K700E. 
with_normalized_insertion <- upstream_by_cb_prop_celltype %>% 
  mutate(normalized_length = upstream_insertion_length * prop) %>% 
  group_by(celltype, cb) %>% 
  summarise(normalized_insertion_length = sum(normalized_length), total_count = mean(total_count))

normalized_insertion_pivot <- with_normalized_insertion %>% filter(total_count > 30) %>% 
  select(-total_count) %>% pivot_wider(names_from = celltype, values_from = normalized_insertion_length, values_fill = 0)

normalized_insertion_pivot_no_zeros <- normalized_insertion_pivot %>% filter(K700E != WT) %>% filter(abs(K700E - WT) > 0.5)
ggplot(normalized_insertion_pivot_no_zeros, aes(WT, K700E)) + geom_point(alpha = 0.3)


### Merge the p-val with the normalized values
with_pval_normalized_insertion <- merge(normalized_insertion_pivot_no_zeros, with_and_without_insertion_merged, by = "cb") %>% 
  as_tibble()

with_pval_normalized_insertion <- merge(with_pval_normalized_insertion, bc_table_all, by.x = "cb", by.y = "barcodeRevcomp") %>%
  separate(ID, sep = ";", into = c("ENSG", "gene", "locus"))

highlight_data2 <- with_pval_normalized_insertion %>% filter(log_p_val > 10) 
ggplot(with_pval_normalized_insertion, aes(K700E, log_p_val)) + geom_point(color = "gray") +  
  geom_point(data = highlight_data2, shape = 21, size = 1.75, color = "#e64c35", fill = "#e64c35") +
  geom_text_repel(data = highlight_data2, aes(label = gene)) + 
  theme_bw()+ scale_fill_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  xlab("K700E Normalized Retained Intron Length") + 
  ylab("-log(p value)")
ggsave(file="~/ManuscriptFigures/22_K700E_pval.pdf", units = "cm", width = 13, height = 12, dpi = 300) #saves g
ggsave(file="~/ManuscriptFigures/22_K700E_pval.png", units = "cm", width = 13, height = 12, dpi = 300) #saves g


ggplot(with_pval_normalized_insertion, aes(WT, K700E)) + geom_point(alpha = 0.5, color = "gray") + 
  geom_point(data = highlight_data2, shape = 21, size = 1.75, color = "#e64c35", fill = "#e64c35") +
  geom_text_repel(data = highlight_data2, aes(label = gene)) +
  theme_bw()+ scale_fill_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave(file="~/ManuscriptFigures/23_K700E_vs_WT.pdf", units = "cm", width = 13, height = 12, dpi = 300) #saves g
ggsave(file="~/ManuscriptFigures/23_K700E_vs_WT.png", units = "cm", width = 13, height = 12, dpi = 300) #saves g
