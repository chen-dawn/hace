library(tidyverse)
library(vroom)
library(data.table)
library(ggpubr)
library(gridExtra)
options(scipen=0)
samples <- c("AID_S11",
             "AID_S12",
             "BLM_S1",
             "BLM_S2",
             "nCas9_S10",
             "nCas9_S9",
             "Ns3h_S3",
             "Ns3h_S4",
             "PcrA-M6_S7",
             "PcrA-M6_S8",
             "PcrA_S5",
             "PcrA_S6",
             "pUC19_S15",
             "pUC19_S16")
out_dir <- "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/"

#### THIS IS 1M CHUNK SIZE #####
all_filepath <- list.files(file.path("/Volumes/broad_thechenlab/Dawn/Helicase/230905ExomeSeq/chunk_stats/"), pattern="merged_1M_stats", full.names = T)
df <- vroom(all_filepath, id= "filename")
df_cleaned <- df %>% mutate(sample = basename(filename)) %>%
  mutate(sample = gsub("merged_1M_stats_","", sample)) %>%
  mutate(sample = gsub(".tsv","", sample)) %>% 
  mutate(condition = str_extract(sample, ".*(?=_S\\d)")) %>% 
  select(-filename)
  
df_cleaned_filtered <- df_cleaned %>% filter(num_base > 0) %>% filter(num_C > 100)

# CtoT by counts for fisher test.
grouped_by_counts <- df_cleaned_filtered %>% group_by(chr, start, end, condition) %>% 
  summarise(num_CtoT_Csum = mean(num_CtoT_Csum), 
            num_CtoT_Tsum = mean(num_CtoT_Tsum), 
            num_CtoT_Asum = mean(num_CtoT_Asum), 
            num_CtoT_Gsum = mean(num_CtoT_Gsum)) %>% 
  mutate(num_CtoT_rate = num_CtoT_Tsum/(num_CtoT_Tsum + num_CtoT_Csum))

CtoT_pivot_by_counts <- grouped_by_counts %>% select(-num_CtoT_Asum, -num_CtoT_Gsum) %>% 
  pivot_wider(names_from = "condition", values_from = starts_with("num"))


CtoT_pivot_by_counts$fisher_res_AID <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_AID']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_AID']]
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

CtoT_pivot_by_counts$fisher_res_BLM <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_BLM']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_BLM']]
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
CtoT_pivot_by_counts$fisher_res_Ns3h <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_Ns3h']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_Ns3h']]
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
CtoT_pivot_by_counts$fisher_res_PcrA <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA']]
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
CtoT_pivot_by_counts$`fisher_res_PcrA-M6` <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA-M6']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA-M6']]
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
  mutate(fisher_res_AID.adj = p.adjust(fisher_res_AID, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_BLM.adj = p.adjust(fisher_res_BLM, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_Ns3h.adj = p.adjust(fisher_res_Ns3h, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_PcrA.adj = p.adjust(fisher_res_PcrA, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(`fisher_res_PcrA-M6.adj` = p.adjust(`fisher_res_PcrA-M6`, method = "BH", n = nrow(CtoT_pivot_by_counts))) 


pAID <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_AID, color = fisher_res_AID.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + scale_color_viridis_c(option = "D")
pBLM <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_BLM, color = fisher_res_BLM.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pNs3h <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_Ns3h, color = fisher_res_Ns3h.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pPcrA <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_PcrA, color = fisher_res_PcrA.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pPcrAM6 <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, `num_CtoT_rate_PcrA-M6`, color = `fisher_res_PcrA-M6.adj`)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, main=text_grob("CtoT By Counts 1M bp chunks"))
as_ggplot(g)

paste("Bin size 1M bp")
high_AID <- CtoT_pivot_by_counts %>% filter(fisher_res_AID.adj < 0.05)
paste("High AID:", nrow(high_AID))

high_BLM <- CtoT_pivot_by_counts %>% filter(fisher_res_BLM.adj < 0.05)
paste("High BLM:", nrow(high_BLM))

high_Ns3h <- CtoT_pivot_by_counts %>% filter(fisher_res_Ns3h.adj < 0.05)
paste("High Ns3h:", nrow(high_Ns3h))

high_PcrA <- CtoT_pivot_by_counts %>% filter(fisher_res_PcrA.adj < 0.05)
paste("High PcrA:", nrow(high_PcrA))

high_PcrAM6 <- CtoT_pivot_by_counts %>% filter(`fisher_res_PcrA-M6.adj` < 0.05)
paste("High PcrAM6:", nrow(high_PcrAM6))
# ggplot(df_cleaned_filtered, aes(CtoA_by_count, CtoA_by_avgMutRate)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~condition)

# CtoT by counts.
CtoTby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_count", values_fn = mean)

pAID <- ggplot(CtoTby_count, aes(nCas9, AID))+ geom_point() + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_count, aes(nCas9, BLM))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_count, aes(nCas9, PcrA))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_count, aes(nCas9, `PcrA-M6`))+ geom_point() + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_count, aes(nCas9, Ns3h))+ geom_point() + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By Counts 1M chunks",gp=gpar(fontsize=20)))
as_ggplot(g)


# CtoT by average mut rate
CtoTby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoTby_avgMutRate, aes(nCas9, AID))+ geom_point() + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_avgMutRate, aes(nCas9, BLM))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_avgMutRate, aes(nCas9, PcrA))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point() + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_avgMutRate, aes(nCas9, Ns3h))+ geom_point() + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By AvgMutRate 1M chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by counts.
CtoAby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_count", values_fn = mean)

pAID <- ggplot(CtoAby_count, aes(nCas9, AID))+ geom_point() + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoAby_count, aes(nCas9, BLM))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoAby_count, aes(nCas9, PcrA))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoAby_count, aes(nCas9, `PcrA-M6`))+ geom_point() + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoAby_count, aes(nCas9, Ns3h))+ geom_point() + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By Counts 1M chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by average mut rate.
CtoAby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoAby_avgMutRate, aes(nCas9, AID))+ geom_point() + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoAby_avgMutRate, aes(nCas9, BLM))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoAby_avgMutRate, aes(nCas9, PcrA))+ geom_point() + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoAby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point() + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoAby_avgMutRate, aes(nCas9, Ns3h))+ geom_point() + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By AvgMutRate 1M chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

##### Try again with 100k chunks #####
all_filepath <- list.files(file.path("/Volumes/broad_thechenlab/Dawn/Helicase/230905ExomeSeq/chunk_stats/"), pattern="merged_100k_stats", full.names = T)
df <- vroom(all_filepath, id= "filename")
df_cleaned <- df %>% mutate(sample = basename(filename)) %>%
  mutate(sample = gsub("merged_100k_stats_","", sample)) %>%
  mutate(sample = gsub(".tsv","", sample)) %>% 
  mutate(condition = str_extract(sample, ".*(?=_S\\d)")) %>% 
  select(-filename)

df_cleaned_filtered <- df_cleaned %>% filter(num_base > 100) %>% filter(num_C > 100) %>% 
  mutate(num_T = num_C * CtoT_by_count/(1-CtoT_by_count)) %>% 
  mutate(num_T = as.integer(num_T))

# CtoT by counts for fisher test.
grouped_by_counts <- df_cleaned_filtered %>% group_by(chr, start, end, condition) %>% 
  summarise(num_CtoT_Csum = mean(num_CtoT_Csum), 
            num_CtoT_Tsum = mean(num_CtoT_Tsum), 
            num_CtoT_Asum = mean(num_CtoT_Asum), 
            num_CtoT_Gsum = mean(num_CtoT_Gsum)) %>% 
  mutate(num_CtoT_rate = num_CtoT_Tsum/(num_CtoT_Tsum + num_CtoT_Csum))

CtoT_pivot_by_counts <- grouped_by_counts %>% select(-num_CtoT_Asum, -num_CtoT_Gsum) %>% 
  pivot_wider(names_from = "condition", values_from = starts_with("num"))


CtoT_pivot_by_counts$fisher_res_AID <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_AID']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_AID']]
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

CtoT_pivot_by_counts$fisher_res_BLM <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_BLM']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_BLM']]
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
CtoT_pivot_by_counts$fisher_res_Ns3h <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_Ns3h']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_Ns3h']]
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
CtoT_pivot_by_counts$fisher_res_PcrA <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA']]
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
CtoT_pivot_by_counts$`fisher_res_PcrA-M6` <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA-M6']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA-M6']]
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
CtoT_pivot_by_counts$fisher_res_pUC19 <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_pUC19']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_pUC19']]
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
  mutate(fisher_res_AID.adj = p.adjust(fisher_res_AID, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_BLM.adj = p.adjust(fisher_res_BLM, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_Ns3h.adj = p.adjust(fisher_res_Ns3h, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_PcrA.adj = p.adjust(fisher_res_PcrA, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(`fisher_res_PcrA-M6.adj` = p.adjust(`fisher_res_PcrA-M6`, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_pUC19.adj = p.adjust(fisher_res_pUC19, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(chr_pos = paste(chr, start, sep="_"))
  

# We filter out some of the regions that are also consistently high in control. 
high_pUC19 <- CtoT_pivot_by_counts %>% filter(fisher_res_pUC19.adj < 0.05)
chr_pos_to_filter <- c(high_pUC19$chr_pos, "chr15_90800001")
CtoT_pivot_by_counts <- CtoT_pivot_by_counts %>% filter(!(chr_pos%in% chr_pos_to_filter))


pAID <- ggplot(CtoT_pivot_by_counts %>% arrange(desc(fisher_res_AID.adj)), aes(log10(num_CtoT_rate_nCas9), log10(num_CtoT_rate_AID), color = log10(fisher_res_AID.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("AID + nCas9") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0), name = "log10(p-value)") 

pBLM <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(fisher_res_BLM.adj)), aes(log10(num_CtoT_rate_nCas9), log10(num_CtoT_rate_BLM), color = log10(fisher_res_BLM.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("AID-BLM-UGI + nCas9") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0), name = "log10(p-value)") 

pNs3h <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(fisher_res_Ns3h.adj)), aes(log10(num_CtoT_rate_nCas9), log10(num_CtoT_rate_Ns3h), color = log10(fisher_res_Ns3h.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("AID-Ns3h-UGI + nCas9") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674",  limits = c(-13, 0),name = "log10(p-value)") 

pPcrA <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(fisher_res_PcrA.adj)), aes(log10(num_CtoT_rate_nCas9), log10(num_CtoT_rate_PcrA), color = log10(fisher_res_PcrA.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("AID-PcrA-UGI + nCas9") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0),name = "log10(p-value)") 

pPcrAM6 <- ggplot(CtoT_pivot_by_counts%>% arrange(desc(`fisher_res_PcrA-M6.adj`)), aes(log10(num_CtoT_rate_nCas9), log10(`num_CtoT_rate_PcrA-M6`), color = log10(`fisher_res_PcrA-M6.adj`))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("AID-PcrA M6-UGI + nCas9") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674",  limits = c(-13, 0),name = "log10(p-value)") 

ppUC19 <- ggplot(CtoT_pivot_by_counts, aes(log10(num_CtoT_rate_nCas9), log10(num_CtoT_rate_pUC19), color = log10(fisher_res_pUC19.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(-5,-2) + ylim(-5,-2) +
  xlab("nCas9 only") + ylab("WT") + 
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", limits = c(-13, 0),name = "log10(p-value)") 

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, ppUC19, nrow = 2)

# as_ggplot(g)
ggsave(file.path(out_dir, "supp3_exome_analysis_redo.pdf"), g, units = "cm", width = 30, height = 13, dpi = 300)

paste("Bin size 100k bp bp")
high_AID <- CtoT_pivot_by_counts %>% filter(fisher_res_AID.adj < 0.05) %>% mutate(high_sample = "AID")
paste("High AID:", nrow(high_AID))

high_BLM <- CtoT_pivot_by_counts %>% filter(fisher_res_BLM.adj < 0.05)  %>% mutate(high_sample = "BLM")
paste("High BLM:", nrow(high_BLM))

high_Ns3h <- CtoT_pivot_by_counts %>% filter(fisher_res_Ns3h.adj < 0.05) %>% mutate(high_sample = "Ns3h")
paste("High Ns3h:", nrow(high_Ns3h))

high_PcrA <- CtoT_pivot_by_counts %>% filter(fisher_res_PcrA.adj < 0.05) %>% mutate(high_sample = "PcrA")
paste("High PcrA:", nrow(high_PcrA))

high_PcrAM6 <- CtoT_pivot_by_counts %>% filter(`fisher_res_PcrA-M6.adj` < 0.05) %>% mutate(high_sample = "PcrA M6")
paste("High PcrAM6:", nrow(high_PcrAM6))

high_pUC19 <- CtoT_pivot_by_counts %>% filter(fisher_res_pUC19.adj < 0.05) %>% mutate(high_sample = "pUC19")
paste("High pUC19:", nrow(high_pUC19))

# look at high in all:
high_all <- rbind(high_AID, high_BLM, high_Ns3h, high_PcrA, high_PcrAM6, high_pUC19) %>% arrange(chr, start) %>% select(chr, start, end, high_sample)

high_grouped <- high_all %>% group_by(chr, start, end) %>% 
  summarise(n = n(), values = paste(high_sample, collapse = ",")) %>% 
  arrange(chr, start)

write.csv(high_grouped, file.path(out_dir, "supp3_exome_analysis_datatable.csv"))



ggplot(CtoT_pivot_by_counts %>% arrange(desc(fisher_res_AID.adj)), aes(log10(num_CtoT_rate_AID/num_CtoT_rate_nCas9), -log10(fisher_res_AID.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", midpoint = 0.8, name = "P-value") 
ggplot(CtoT_pivot_by_counts %>% arrange(desc(fisher_res_BLM.adj)), aes(log10(num_CtoT_rate_BLM/num_CtoT_rate_nCas9), -log10(fisher_res_BLM.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", midpoint = 0.8, name = "P-value") 
ggplot(CtoT_pivot_by_counts %>% arrange(desc(fisher_res_Ns3h.adj)), aes(log10(num_CtoT_rate_Ns3h/num_CtoT_rate_nCas9), -log10(fisher_res_Ns3h.adj))) + 
  geom_point(alpha = 0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_gradient2(high = "#E9F6E5",mid = "#5DD0BE",low = "#003674", midpoint = 0.8, name = "P-value") 

# ggplot(df_cleaned_filtered, aes(CtoT_by_count, CtoT_by_avgMutRate)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~condition)

# CtoT by counts.
CtoTby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_count", values_fn = mean)

pAID <- ggplot(CtoTby_count, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_count, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_count, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_count, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_count, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By Counts 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoT by average mut rate
CtoTby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoTby_avgMutRate, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_avgMutRate, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_avgMutRate, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_avgMutRate, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By AvgMutRate 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by counts.
CtoAby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_count", values_fn = mean)

pAID <- ggplot(CtoAby_count, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoAby_count, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoAby_count, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoAby_count, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoAby_count, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By Counts 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by average mut rate.
CtoAby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoAby_avgMutRate, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10() + geom_abline(slope=1, intercept = 0)
pBLM <- ggplot(CtoAby_avgMutRate, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrA <- ggplot(CtoAby_avgMutRate, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrAM6 <- ggplot(CtoAby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pNs3h <- ggplot(CtoAby_avgMutRate, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By AvgMutRate 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)


CtoAby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, sample, CtoA_by_avgMutRate) %>% 
  pivot_wider(names_from = "sample", values_from = "CtoA_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoAby_avgMutRate, aes(AID_S11, AID_S12))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10() + geom_abline(slope=1, intercept = 0)
pBLM <- ggplot(CtoAby_avgMutRate, aes(BLM_S1, BLM_S2))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrA <- ggplot(CtoAby_avgMutRate, aes(PcrA_S5, PcrA_S6))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrAM6 <- ggplot(CtoAby_avgMutRate, aes(`PcrA-M6_S7`, `PcrA-M6_S8`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pNs3h <- ggplot(CtoAby_avgMutRate, aes(Ns3h_S3, Ns3h_S4))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By AvgMutRate 100k chunks reps",gp=gpar(fontsize=20)))
as_ggplot(g)


##### Try again with 10k chunks #####
all_filepath <- list.files(file.path("/Volumes/broad_thechenlab/Dawn/Helicase/230905ExomeSeq/chunk_stats/"), pattern="merged_10k_stats", full.names = T)
df <- vroom(all_filepath, id= "filename")
df_cleaned <- df %>% mutate(sample = basename(filename)) %>%
  mutate(sample = gsub("merged_10k_stats_","", sample)) %>%
  mutate(sample = gsub(".tsv","", sample)) %>% 
  mutate(condition = str_extract(sample, ".*(?=_S\\d)")) %>% 
  select(-filename)

df_cleaned_filtered <- df_cleaned %>% filter(num_base > 100) %>% filter(num_C > 50) 

# CtoT by counts for fisher test.
grouped_by_counts <- df_cleaned_filtered %>% group_by(chr, start, end, condition) %>% 
  summarise(num_CtoT_Csum = mean(num_CtoT_Csum), 
            num_CtoT_Tsum = mean(num_CtoT_Tsum), 
            num_CtoT_Asum = mean(num_CtoT_Asum), 
            num_CtoT_Gsum = mean(num_CtoT_Gsum)) %>% 
  mutate(num_CtoT_rate = num_CtoT_Tsum/(num_CtoT_Tsum + num_CtoT_Csum))

CtoT_pivot_by_counts <- grouped_by_counts %>% select(-num_CtoT_Asum, -num_CtoT_Gsum) %>% 
  pivot_wider(names_from = "condition", values_from = starts_with("num"))


CtoT_pivot_by_counts$fisher_res_AID <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_AID']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_AID']]
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

CtoT_pivot_by_counts$fisher_res_BLM <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_BLM']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_BLM']]
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
CtoT_pivot_by_counts$fisher_res_Ns3h <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_Ns3h']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_Ns3h']]
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
CtoT_pivot_by_counts$fisher_res_PcrA <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA']]
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
CtoT_pivot_by_counts$`fisher_res_PcrA-M6` <- apply(CtoT_pivot_by_counts, 1, function(x){
  p1 <- x[['num_CtoT_Csum_nCas9']]
  p2 <- x[['num_CtoT_Csum_PcrA-M6']]
  p3 <- x[['num_CtoT_Tsum_nCas9']]
  p4 <- x[['num_CtoT_Tsum_PcrA-M6']]
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
  mutate(fisher_res_AID.adj = p.adjust(fisher_res_AID, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_BLM.adj = p.adjust(fisher_res_BLM, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_Ns3h.adj = p.adjust(fisher_res_Ns3h, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(fisher_res_PcrA.adj = p.adjust(fisher_res_PcrA, method = "BH", n = nrow(CtoT_pivot_by_counts))) %>%
  mutate(`fisher_res_PcrA-M6.adj` = p.adjust(`fisher_res_PcrA-M6`, method = "BH", n = nrow(CtoT_pivot_by_counts))) 


pAID <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_AID, color = fisher_res_AID.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10() + scale_color_viridis_c(option = "D")
pBLM <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_BLM, color = fisher_res_BLM.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pNs3h <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_Ns3h, color = fisher_res_Ns3h.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pPcrA <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, num_CtoT_rate_PcrA, color = fisher_res_PcrA.adj)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
pPcrAM6 <- ggplot(CtoT_pivot_by_counts, aes(num_CtoT_rate_nCas9, `num_CtoT_rate_PcrA-M6`, color = `fisher_res_PcrA-M6.adj`)) + geom_point(alpha = 0.2) + scale_x_log10() + scale_y_log10()+ scale_color_viridis_c(option = "D")
g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, main=text_grob("CtoT By Counts 10k chunks"))
as_ggplot(g)

paste("Bin size 10k bp")
high_AID <- CtoT_pivot_by_counts %>% filter(fisher_res_AID.adj < 0.05)
paste("High AID:", nrow(high_AID))

high_BLM <- CtoT_pivot_by_counts %>% filter(fisher_res_BLM.adj < 0.05)
paste("High BLM:", nrow(high_BLM))

high_Ns3h <- CtoT_pivot_by_counts %>% filter(fisher_res_Ns3h.adj < 0.05)
paste("High Ns3h:", nrow(high_Ns3h))

high_PcrA <- CtoT_pivot_by_counts %>% filter(fisher_res_PcrA.adj < 0.05)
paste("High PcrA:", nrow(high_PcrA))

high_PcrAM6 <- CtoT_pivot_by_counts %>% filter(`fisher_res_PcrA-M6.adj` < 0.05)
paste("High PcrAM6:", nrow(high_PcrAM6))

# ggplot(df_cleaned_filtered, aes(CtoT_by_count, CtoT_by_avgMutRate)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~condition)

# CtoT by counts.
CtoTby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_count", values_fn = mean)

pAID <- ggplot(CtoTby_count, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_count, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_count, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_count, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_count, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By Counts 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoT by average mut rate
CtoTby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoT_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoT_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoTby_avgMutRate, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoTby_avgMutRate, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoTby_avgMutRate, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoTby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoTby_avgMutRate, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoT By AvgMutRate 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by counts.
CtoAby_count <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_count) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_count", values_fn = mean)

pAID <- ggplot(CtoAby_count, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pBLM <- ggplot(CtoAby_count, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrA <- ggplot(CtoAby_count, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pPcrAM6 <- ggplot(CtoAby_count, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()
pNs3h <- ggplot(CtoAby_count, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By Counts 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)

# CtoA by average mut rate.
CtoAby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, condition, CtoA_by_avgMutRate) %>% 
  pivot_wider(names_from = "condition", values_from = "CtoA_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoAby_avgMutRate, aes(nCas9, AID))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10() + geom_abline(slope=1, intercept = 0)
pBLM <- ggplot(CtoAby_avgMutRate, aes(nCas9, BLM))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrA <- ggplot(CtoAby_avgMutRate, aes(nCas9, PcrA))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrAM6 <- ggplot(CtoAby_avgMutRate, aes(nCas9, `PcrA-M6`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pNs3h <- ggplot(CtoAby_avgMutRate, aes(nCas9, Ns3h))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By AvgMutRate 100k chunks",gp=gpar(fontsize=20)))
as_ggplot(g)


CtoAby_avgMutRate <- df_cleaned_filtered %>% select(chr, start, end, sample, CtoA_by_avgMutRate) %>% 
  pivot_wider(names_from = "sample", values_from = "CtoA_by_avgMutRate", values_fn = mean)

pAID <- ggplot(CtoAby_avgMutRate, aes(AID_S11, AID_S12))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10() + geom_abline(slope=1, intercept = 0)
pBLM <- ggplot(CtoAby_avgMutRate, aes(BLM_S1, BLM_S2))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrA <- ggplot(CtoAby_avgMutRate, aes(PcrA_S5, PcrA_S6))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pPcrAM6 <- ggplot(CtoAby_avgMutRate, aes(`PcrA-M6_S7`, `PcrA-M6_S8`))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)
pNs3h <- ggplot(CtoAby_avgMutRate, aes(Ns3h_S3, Ns3h_S4))+ geom_point(color = "salmon") + scale_x_log10() + scale_y_log10()+ geom_abline(slope=1, intercept = 0)

g <- arrangeGrob(pAID, pBLM, pNs3h, pPcrA, pPcrAM6, top=textGrob("CtoA By AvgMutRate 100k chunks reps",gp=gpar(fontsize=20)))
as_ggplot(g)
