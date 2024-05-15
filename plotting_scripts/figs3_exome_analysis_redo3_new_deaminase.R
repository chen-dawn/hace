library(tidyverse)
library(vroom)
library(data.table)
library(ggpubr)
library(gridExtra)
options(scipen=0)

# Helper function for Fisher's Exact Test
fisher_test_helper <- function(p1, p2, p3, p4) {
  if (is.na(p1) | is.na(p2) | is.na(p3) | is.na(p4)) {
    return(NA)
  }
  
  dat <- matrix(c(as.integer(p1), as.integer(p3), as.integer(p2), as.integer(p4)), nrow = 2)
  test <- fisher.test(dat)
  return(test$p.value)
}

samples <- c("Exome-DC391_S1","Exome-DC391_S5","Exome-DC646_S2","Exome-DC646_S6","Exome-pBO127_S3","Exome-pBO127_S7","Exome-pUC19_S4","Exome-pUC19_S8")
out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"

##### Try again with 100k chunks #####
all_filepath <- list.files(file.path("/Volumes/broad_thechenlab/Dawn/Helicase/240417ExomeSeq/combined_stats_100k/"), pattern="merged_100k_stats", full.names = T)
df <- vroom(all_filepath, id= "filename")
df_cleaned <- df %>% mutate(sample = basename(filename)) %>%
  mutate(sample = gsub("merged_100k_stats_","", sample)) %>%
  mutate(sample = gsub(".tsv","", sample)) %>% 
  mutate(sample=gsub("Exome-", "", sample)) %>% 
  mutate(condition = str_extract(sample, ".*(?=_S\\d)")) %>% 
  mutate(condition = gsub("Exome-", "", condition)) %>%
  select(-filename)

df_cleaned_filtered <- df_cleaned %>% 
  mutate(num_base = as.integer(num_base)) %>%
  filter(num_base > 100) %>% 
  filter(num_C > 50) %>% 
  mutate(num_T = num_C * CtoT_by_count/(1-CtoT_by_count)) %>% 
  mutate(num_T = as.integer(num_T)) %>% 
  # I'm going to group the CtoT Csum and GtoA Gsum by counts
  mutate(num_CtoT_Csum = num_CtoT_Csum + num_GtoA_Gsum) %>%
  mutate(num_CtoT_Tsum = num_CtoT_Tsum + num_GtoA_Asum) 

# CtoT by counts for fisher test.
grouped_by_counts <- df_cleaned_filtered %>% group_by(chr, start, end, condition) %>% 
  summarise(num_CtoT_Csum = mean(num_CtoT_Csum), 
            num_CtoT_Tsum = mean(num_CtoT_Tsum), 
            num_CtoT_Asum = mean(num_CtoT_Asum), 
            num_CtoT_Gsum = mean(num_CtoT_Gsum)) %>% 
  mutate(num_CtoT_rate = num_CtoT_Tsum/(num_CtoT_Tsum + num_CtoT_Csum))

CtoT_pivot_by_counts <- grouped_by_counts %>% select(-num_CtoT_Asum, -num_CtoT_Gsum) %>% 
  pivot_wider(names_from = "condition", values_from = starts_with("num"))


CtoT_pivot_by_counts$fisher_res_DC391 <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_pUC19']]
  p2 <- x[['num_CtoT_Csum_DC391']]
  p3 <- x[['num_CtoT_Tsum_pUC19']]
  p4 <- x[['num_CtoT_Tsum_DC391']]
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

CtoT_pivot_by_counts$fisher_res_DC646 <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_pUC19']]
  p2 <- x[['num_CtoT_Csum_DC646']]
  p3 <- x[['num_CtoT_Tsum_pUC19']]
  p4 <- x[['num_CtoT_Tsum_DC646']]
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
CtoT_pivot_by_counts$fisher_res_pBO127 <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_pUC19']]
  p2 <- x[['num_CtoT_Csum_pBO127']]
  p3 <- x[['num_CtoT_Tsum_pUC19']]
  p4 <- x[['num_CtoT_Tsum_pBO127']]
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

CtoT_pivot_by_counts <- CtoT_pivot_by_counts %>% 
  mutate(fisher_res_DC391.adj = p.adjust(fisher_res_DC391, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_DC646.adj = p.adjust(fisher_res_DC646, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_pBO127.adj = p.adjust(fisher_res_pBO127, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(chr_pos = paste(chr, start, sep="_"))
  

# We filter out some of the regions that are also consistently high in control. 
# high_pUC19 <- CtoT_pivot_by_counts %>% filter(fisher_res_pUC19.adj < 0.05)
# chr_pos_to_filter <- c(high_pUC19$chr_pos, "chr15_90800001")
# CtoT_pivot_by_counts <- CtoT_pivot_by_counts %>% filter(!(chr_pos%in% chr_pos_to_filter))


pDC391 <- ggplot(CtoT_pivot_by_counts %>% arrange(desc(fisher_res_DC391.adj)), aes(log10(num_CtoT_rate_pUC19), log10(num_CtoT_rate_DC391), color = log10(fisher_res_DC391.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("pUC19") + ylab("DC391") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0), name = "log10(p-value)") 

pDC646 <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(fisher_res_DC646.adj)), aes(log10(num_CtoT_rate_pUC19), log10(num_CtoT_rate_DC646), color = log10(fisher_res_DC646.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("pUC19") + ylab("DC646") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0), name = "log10(p-value)") 

pBO127 <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(fisher_res_pBO127.adj)), aes(log10(num_CtoT_rate_pUC19), log10(num_CtoT_rate_pBO127), color = log10(fisher_res_pBO127.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("pUC19") + ylab("pBO127") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674",  limits = c(-13, 0),name = "log10(p-value)") 

g <- arrangeGrob(pBO127, pDC646, pDC391, nrow = 1)

# as_ggplot(g)
ggsave(file.path(out_dir, "supp3_CtoT_exome_analysis_redo_revision_2.pdf"), g, units = "cm", width = 30, height = 6.5, dpi = 300)

paste("Bin size 100k bp bp")
high_DC391 <- CtoT_pivot_by_counts %>% filter(fisher_res_DC391.adj < 0.05) %>% mutate(high_sample = "DC391")
paste("High DC391:", nrow(high_DC391))

high_DC646 <- CtoT_pivot_by_counts %>% filter(fisher_res_DC646.adj < 0.05)  %>% mutate(high_sample = "DC646")
paste("High DC646:", nrow(high_DC646))

high_pBO127 <- CtoT_pivot_by_counts %>% filter(fisher_res_pBO127.adj < 0.05) %>% mutate(high_sample = "pBO127")
paste("High pBO127:", nrow(high_pBO127))

# look at high in all:
high_all <- rbind(high_DC391, high_DC646, high_pBO127) %>% arrange(chr, start) %>% select(chr, start, end, high_sample)

high_grouped <- high_all %>% group_by(chr, start, end) %>% 
  summarise(n = n(), values = paste(high_sample, collapse = ",")) %>% 
  arrange(chr, start)

write.csv(high_grouped, file.path(out_dir, "supp3_CtoT_exome_analysis_datatable_revision.csv"))


#############
## AtoG
#############
df_cleaned_filtered <- df_cleaned %>% 
  mutate(num_base = as.integer(num_base)) %>%
  filter(num_base > 100) %>% 
  filter(num_C > 50) %>% 
  # I'm going to group the CtoT Csum and GtoA Gsum by counts
  mutate(num_AtoG_Asum = num_AtoG_Asum + num_TtoC_Tsum) %>%
  mutate(num_AtoG_Gsum = num_AtoG_Gsum + num_TtoC_Csum) 

grouped_by_counts <- df_cleaned_filtered %>% group_by(chr, start, end, condition) %>% 
  summarise(num_AtoG_Asum = mean(num_AtoG_Asum), 
            num_AtoG_Gsum = mean(num_AtoG_Gsum)) %>% 
  mutate(num_AtoG_rate = num_AtoG_Gsum/(num_AtoG_Gsum + num_AtoG_Asum))

AtoG_pivot_by_counts <- grouped_by_counts %>% 
  pivot_wider(names_from = "condition", values_from = starts_with("num"))

AtoG_pivot_by_counts$fisher_res_DC391 <- apply(AtoG_pivot_by_counts, 1, function(x){
  p1 <- x[['num_AtoG_Asum_pUC19']]
  p2 <- x[['num_AtoG_Asum_DC391']]
  p3 <- x[['num_AtoG_Gsum_pUC19']]
  p4 <- x[['num_AtoG_Gsum_DC391']]
  if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
    return(NA)
  }
  
  dat <- matrix(c(as.integer(p1), as.integer(p3), as.integer(p2), as.integer(p4)), nrow = 2)
  test <- fisher.test(dat)
  return(test$p.value)
})

AtoG_pivot_by_counts$fisher_res_DC646 <- apply(AtoG_pivot_by_counts, 1, function(x){
  p1 <- x[['num_AtoG_Asum_pUC19']]
  p2 <- x[['num_AtoG_Asum_DC646']]
  p3 <- x[['num_AtoG_Gsum_pUC19']]
  p4 <- x[['num_AtoG_Gsum_DC646']]
  if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
    return(NA)
  }
  
  dat <- matrix(c(as.integer(p1), as.integer(p3), as.integer(p2), as.integer(p4)), nrow = 2)
  test <- fisher.test(dat)
  return(test$p.value)
})

AtoG_pivot_by_counts$fisher_res_pBO127 <- apply(AtoG_pivot_by_counts, 1, function(x){
  p1 <- x[['num_AtoG_Asum_pUC19']]
  p2 <- x[['num_AtoG_Asum_pBO127']]
  p3 <- x[['num_AtoG_Gsum_pUC19']]
  p4 <- x[['num_AtoG_Gsum_pBO127']]
  if(is.na(p1)|is.na(p2)|is.na(p3)|is.na(p4)){
    return(NA)
  }
  
  dat <- matrix(c(as.integer(p1), as.integer(p3), as.integer(p2), as.integer(p4)), nrow = 2)
  test <- fisher.test(dat)
  return(test$p.value)
})

AtoG_pivot_by_counts <- AtoG_pivot_by_counts %>% 
  mutate(fisher_res_DC391.adj = p.adjust(fisher_res_DC391, method = "BH", n = nrow(AtoG_pivot_by_counts))) %>%
  mutate(fisher_res_DC646.adj = p.adjust(fisher_res_DC646, method = "BH", n = nrow(AtoG_pivot_by_counts))) %>%
  mutate(fisher_res_pBO127.adj = p.adjust(fisher_res_pBO127, method = "BH", n = nrow(AtoG_pivot_by_counts))) %>%
  mutate(chr_pos = paste(chr, start, sep="_"))


# Visualization of A to G transition rates
pDC391 <- ggplot(AtoG_pivot_by_counts %>% arrange(desc(fisher_res_DC391.adj)), aes(log10(num_AtoG_rate_pUC19), log10(num_AtoG_rate_DC391), color = log10(fisher_res_DC391.adj))) +
  geom_point(alpha = 0.3) + theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  xlim(-5, -2) + ylim(-5, -2) +
  xlab("pUC19") + ylab("DC391") +
  scale_colour_gradient2(high = "#E9F6E5", mid = "#5DD0BE", low = "#003674", limits = c(-13, 0), name = "log10(p-value)")

pDC646 <- ggplot(AtoG_pivot_by_counts %>% arrange(desc(fisher_res_DC646.adj)), aes(log10(num_AtoG_rate_pUC19), log10(num_AtoG_rate_DC646), color = log10(fisher_res_DC646.adj))) +
  geom_point(alpha = 0.3) + theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  xlim(-5, -2) + ylim(-5, -2) +
  xlab("pUC19") + ylab("DC646") +
  scale_colour_gradient2(high = "#E9F6E5", mid = "#5DD0BE", low = "#003674", limits = c(-13, 0), name = "log10(p-value)")

pBO127 <- ggplot(AtoG_pivot_by_counts %>% arrange(desc(fisher_res_pBO127.adj)), aes(log10(num_AtoG_rate_pUC19), log10(num_AtoG_rate_pBO127), color = log10(fisher_res_pBO127.adj))) +
  geom_point(alpha = 0.3) + theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  xlim(-5, -2) + ylim(-5, -2) +
  xlab("pUC19") + ylab("rAPOBEC1") +
  scale_colour_gradient2(high = "#E9F6E5", mid = "#5DD0BE", low = "#003674", limits = c(-13, 0), name = "log10(p-value)")

# Arrange plots in a single row
g <- arrangeGrob(pBO127, pDC646, pDC391, nrow = 1)

# Save plot to file
ggsave(file.path(out_dir, "supp3_AtoG_exome_analysis_revision_2.pdf"), g, units = "cm", width = 30, height = 6.5, dpi = 300)

# Extract and output significant results
high_DC391 <- AtoG_pivot_by_counts %>% filter(fisher_res_DC391.adj < 0.05) %>% mutate(high_sample = "DC391")
high_DC646 <- AtoG_pivot_by_counts %>% filter(fisher_res_DC646.adj < 0.05) %>% mutate(high_sample = "DC646")
high_pBO127 <- AtoG_pivot_by_counts %>% filter(fisher_res_pBO127.adj < 0.05) %>% mutate(high_sample = "pBO127")

# Combine high results from all conditions
high_all <- rbind(high_DC391, high_DC646, high_pBO127) %>% arrange(chr, start) %>% select(chr, start, end, high_sample)

# Group by genomic region and summarize
high_grouped <- high_all %>% group_by(chr, start, end) %>% 
  summarise(n = n(), values = paste(high_sample, collapse = ",")) %>% 
  arrange(chr, start)

# Write the results to a CSV file
write.csv(high_grouped, file.path(out_dir, "supp3_AtoG_exome_analysis_datatable_revision.csv"))

