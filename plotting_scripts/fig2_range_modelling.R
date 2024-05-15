library(zoo)
library(plyr)
library(vroom)
library(tidyverse)
library(data.table)
library(gridExtra)
library(viridis)
library(ggsci)
library(ggpubr)
library(cowplot)
library(ggforce)


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
 
# MEK1_dir <- "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/"
# MEK1_pileup_filenames <- list.files(file.path(MEK1_dir), pattern="pileup.tsv", full.names = T)
# 
# all_files <- vroom(MEK1_pileup_filenames, id = "filename")
# all_files_df <- all_files %>%  mutate(sample = str_extract(basename(filename), ".+(S\\d+)")) %>% 
#   mutate(condition = str_extract(sample, "^.+(?=-[:upper:]\\d+_S\\d+)")) %>% 
#   select(-filename)
# 
# all_files_df <- all_files_df %>% 
#   mutate(condition = ifelse(is.na(condition), str_extract(sample, "^.+(?=_\\d_S\\d+)"), condition)) %>% 
#   mutate(condition = gsub("_", "-", condition))
# 
# fwrite(all_files_df, "~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")

all_files_df <- fread("~/Dropbox (Harvard University)/03Helicase/data/Expt06_Multiple_locations_range_tagment/merged_all_tagment_pileups.csv")
MEK1_pileup <- all_files_df %>% 
  mutate(base_sum = (A + C + `T` + G )) 
#  filter(base_sum > 10000)
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
DC203 <- MEK1_pileup %>% filter(grepl("DC203", sample)) %>% 
  mutate(n_base = n_base - 48027) %>% filter(chr == "10_MAP2K1") %>% mutate(loci = "MAP2K1")

RUNX1 <- MEK1_pileup %>% filter(grepl("RUNX1", sample)) %>% 
  mutate(n_base = n_base - 5115) %>% filter(chr == "RUNX1_8kb_extraction") %>% mutate(loci = "RUNX1")
RUNX1 <- flip_DNA(RUNX1)

DC329 <- MEK1_pileup %>% filter(grepl("DC329", sample)) %>% 
  mutate(n_base = n_base - 2478) %>% filter(chr == "chr6_(NC_000006)_TNF_long_extraction") %>% mutate(loci = "TNF")
DC329 <- flip_DNA(DC329)

DC330 <- MEK1_pileup %>% filter(grepl("DC330", sample)) %>% 
  mutate(n_base = n_base - 1534) %>% filter(chr == "NC_000007_-_IL6_gene") %>% mutate(loci = "IL6")
DC330 <- flip_DNA(DC330)

DC331 <- MEK1_pileup %>% filter(grepl("DC331", sample)) %>% 
  mutate(n_base = n_base - 15562) %>% filter(chr =="NC_000014_-_AKT1_gene") %>% mutate(loci = "AKT1")

DC573 <- MEK1_pileup %>% filter(grepl("DC573", sample)) %>% 
  mutate(n_base = n_base - 44231) %>% filter(chr == "DC573_02_HEK3") %>% mutate(loci = "HEK3")

DC576 <- MEK1_pileup %>% filter(grepl("DC576", sample)) %>% 
  mutate(n_base = n_base - 2972) %>% filter(chr == "DC576_05_VEGFA") %>% mutate(loci = "VEGFA")
DC576 <- flip_DNA(DC576)

DC577 <- MEK1_pileup %>% filter(grepl("DC577", sample)) %>% 
  mutate(n_base = n_base - 7846) %>% filter(chr == "DC577_06_CD209") %>% mutate(loci = "CD209")
DC577 <- flip_DNA(DC577)

# Need to change names.
reconcat_all <- rbind(RUNX1, DC203, DC329, DC330, DC331, DC573, DC576, DC577) %>% filter(n_base > -2000 & n_base < 2000)
# with_names_edited <- reconcat_all %>% mutate(temp_name = str_extract(sample, "^.*(?=_\\d_S\\d+)")) %>%
#   mutate(temp_extension = ifelse(is.na(temp_name),NA, str_extract(sample, "_\\d_S\\d+"))) %>% 
#   mutate(temp_name = gsub("_", "-",temp_name)) %>% 
#   mutate(condition = ifelse(is.na(temp_name), condition, temp_name)) %>% 
#   select(-temp_name, -temp_extension)

# And then now I just extract the guide, nCas9, and helicase from each sample. 
with_plasmid_annotation <- reconcat_all %>% 
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
  mutate(isControl = ifelse(condition %in% c("RUNX1-Ctrl", "DC203-Control", "DC329-Ctrl", "DC330-Ctrl", "DC331-Ctrl", "DC573-DC573", "DC576-DC573", "DC577-DC573"), 1, 0))
  
control_rates <- with_mutation_rate %>% filter(isControl  ==  1) %>% 
  group_by(chr, n_base, ref_base)%>%
  dplyr::summarise(meanCtoT_ctrl = mean(CtoT), 
            meanGtoA_ctrl = mean(GtoA))

with_mut_and_control <- merge(with_mutation_rate, control_rates, by = c("chr", "n_base"), all.x = T) %>% as_tibble() %>% arrange(n_base)
top_ctrl <- with_mut_and_control %>% filter(meanGtoA_ctrl < 5)

with_mut_background_subtracted <- with_mut_and_control %>% 
  filter(meanGtoA_ctrl < 5) %>% 
  mutate(CtoT_bg_subtracted = CtoT) %>% 
  mutate(GtoA_bg_subtracted = GtoA) %>%
  mutate(CtoT_bg_subtracted = CtoT - meanCtoT_ctrl) %>% 
  mutate(GtoA_bg_subtracted = GtoA - meanGtoA_ctrl) %>%
  mutate(CtoA = A/(C + A)) %>%
  mutate(TtoA = A/(A + `T`)) %>%
  filter(n_base >= -1500 & n_base <= 1500) %>% 
  filter(loci %in% c("IL6", "TNF", "HEK3")) %>% 
  filter(!is.na(Cas9Display))


avg_mut_per_base <- with_mut_background_subtracted %>% 
  pivot_longer(cols = contains("bg_subtracted"), names_to = "editMode") %>% 
  group_by(chr, n_base, ref_base.x, editMode, condition, loci, helicaseDisplay, Cas9Display) %>% 
  dplyr::summarise(editRate = mean(value), sd = sd(value)) 

# avg_mut_per_base_filtered <- avg_mut_per_base %>% 
#   filter(ref_base.x %in% c("C", "G"))

##### FIRST FOR FORWARD DIRECTION ######
grouped_by_condition <- with_mut_background_subtracted %>% 
  group_by(condition, n_base, chr, ref_base.y, loci, helicaseDisplay, Cas9Display) %>% 
  dplyr::summarise(GtoA = mean(GtoA_bg_subtracted), CtoT = mean(CtoT_bg_subtracted), GtoA_sd = sd(GtoA_bg_subtracted), CtoT_sd = sd(CtoT_bg_subtracted)) %>% 
  filter(Cas9Display %in% c("nCas9 D10A","nCas9 H840A")| is.na(Cas9Display)) %>% 
  filter(n_base >= 7 & n_base < 1000) %>% 
  filter(GtoA <20 & GtoA > 0)

##### Get Rolling Median #####
unique_condition <- unique(grouped_by_condition$condition)
all_bases <- data.frame(n_base = seq(0, 1500))
my_data <- list()
for (i in seq_along(unique_condition)){
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(condition == unique_condition[i]) 
  condition_name = tmp_mat$condition[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(condition = condition_name) %>% 
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
  mutate(logRollingGtoA = log(rollingGtoA)) 

spacing100_nbase <- seq(100,900, by = 100)
forward_points <- data_for_fit_forward %>% filter(n_base %in% spacing100_nbase)
ggplot(forward_points, aes(n_base, log10(rollingGtoA), color = loci)) + geom_point() + 
  geom_line() + 
  facet_grid(helicaseDisplay~Cas9Display) + ggtitle("101 rolling window width, 100 bp is one point")
write_csv(forward_points, "~/Dropbox (Harvard University)/03Helicase/ManuscriptDraft/R_scripts_outputs/data_supp2_rollingGrouped100bp.csv")

forward_subsample_grouped <- forward_points %>% 
  ungroup() %>% 
  group_by(n_base, helicaseDisplay, Cas9Display) %>% 
  ## summarise(num = n())
  dplyr::summarise(meanLogRolling = mean(log10(rollingGtoA)), 
            sdLogRolling = sd(log10(rollingGtoA)))

ggplot(forward_subsample_grouped, aes(n_base, meanLogRolling)) + 
  geom_point() + 
  geom_line() + 
  facet_grid(helicaseDisplay~Cas9Display) + ggtitle("101 rolling window width, 100 bp is one point")
# If the fit is just on the raw data
# data_for_fit <- grouped_by_condition %>% filter(n_base <= 1000 & !is.na(helicaseDisplay)) %>% 
#   filter(loci != "AKT1") %>% 
#   mutate(logRollingGtoA = log(GtoA)) 


# lm_eqn = function(df){
#   m = lm(logRollingGtoA ~ n_base, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(as.numeric(coef(m)[1]), digits = 2), 
#                         b = format(as.numeric(coef(m)[2]), digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# eq <- ddply(data_for_fit_forward,.(helicaseDisplay, Cas9Display),lm_eqn)


# ggplot(all_roll_means %>% filter(n_base <= 1000 & !is.na(helicaseDisplay)), aes(n_base, rollingGtoA, color = loci)) + geom_point() + 
#   facet_wrap(~helicaseDisplay+Cas9Display) 

# Group by nCas9, helicase, and n_base
grouped_by_base <- data_for_fit_forward %>% 
  dplyr::group_by(helicaseDisplay, Cas9Display, n_base) %>% 
  dplyr::summarise(meanLogRollingGtoA = log(mean(rollingGtoA))) %>% 
  dplyr::mutate(log_n_base = log(n_base))
lm_eqn = function(df){
  m = lm(meanLogRollingGtoA ~ n_base, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(as.numeric(coef(m)[1]), digits = 5), 
                        b = format(as.numeric(coef(m)[2]), digits = 5), 
                        r2 = format(summary(m)$r.squared, digits = 5)))
  as.character(as.expression(eq));                 
}
eq <- ddply(grouped_by_base,.(helicaseDisplay, Cas9Display),lm_eqn)

newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399",  "#FBBF24", "#A3E635", "#E879F9", "#A78BFA")
pD10A <- ggplot(grouped_by_base %>% filter(Cas9Display == "nCas9 D10A"), aes(n_base, meanLogRollingGtoA, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") + ylim(-5,0) + theme(legend.position = "none")
pH840A <- ggplot(grouped_by_base %>% filter(Cas9Display == "nCas9 H840A"), aes(n_base, meanLogRollingGtoA, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") + ylim(-5,0)+ theme(legend.position = c(0.9, 0.8)) 

g <- arrangeGrob(pD10A, pH840A, nrow = 1)
ggsave(file.path(out_dir, "fig2e_range_modelling_helicase.pdf"), g, units = "cm", width = 16, height = 6, dpi = 150)


max_first  <- max(data_for_fit_forward$GtoA, na.rm = T)   # Specify max of first y axis
max_second <- max(data_for_fit_forward$logRollingGtoA, na.rm = T) # Specify max of second y axis
min_first  <- min(data_for_fit_forward$GtoA, na.rm = T)   # Specify min of first y axis
min_second <- -5 # min(data_for_fit_forward$logRollingGtoA, na.rm = T) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}


# Plot the rolling average above the points?
ggplot(data_for_fit_forward %>% filter(loci == "HEK3" & Cas9Display == "nCas9 H840A")) + 
  geom_point(aes(n_base, GtoA)) + 
  facet_wrap(~condition, ncol = 1) + 
  geom_line(aes(n_base, inv_scale_function(logRollingGtoA, scale, shift))) + 
  scale_y_continuous(limits = c(min_first, max_first), sec.axis = sec_axis(~scale_function(., scale, shift), name="logRollingGtoA"))
  




g2 <- ggplot(grouped_by_base, aes(n_base, meanLogRollingGtoA, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") + theme(legend.position = "none") + 
  geom_text(data=eq,aes(x = 300, y = -4,label=V1), parse = TRUE, inherit.aes=FALSE, size = 5) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y~ x) + 
  facet_grid(helicaseDisplay~Cas9Display)

ggsave(file.path(out_dir, "fig2e_range_modelling2.pdf"), g2, units = "cm", width = 30, height = 30, dpi = 300)


# ##### NOW FOR REVERSE DIRECTION ######
# grouped_by_condition <- with_mut_background_subtracted %>% 
#   group_by(condition, n_base, chr, ref_base.y, loci, helicaseDisplay, Cas9Display) %>% 
#   dplyr::summarise(GtoA = mean(GtoA_bg_subtracted), CtoT = mean(CtoT_bg_subtracted)) %>% 
#   filter(Cas9Display %in% c("nCas9 D10A","nCas9 H840A")) %>% 
#   filter(n_base <= -3 & n_base > -1000) %>% 
#   filter(GtoA <20 & GtoA > 0)
# 
# ##### Get Rolling Median #####
# unique_condition <- unique(grouped_by_condition$condition)
# all_bases <- data.frame(n_base = seq(-1000, 0))
# my_data <- list()
# for (i in seq_along(unique_condition)){
#   print(paste("Current i", i))
#   tmp_mat <- grouped_by_condition %>% filter(condition == unique_condition[i]) 
#   condition_name = tmp_mat$condition[1]
#   loci_name = tmp_mat$loci[1]
#   helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
#   Cas9Display_name = tmp_mat$Cas9Display[1]
#   
#   tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
#     as_tibble() %>%
#     filter(n_base < -3 | n_base > 17) %>% 
#     mutate(condition = condition_name) %>% 
#     mutate(loci = loci_name) %>% 
#     mutate(helicaseDisplay = helicaseDisplay_name) %>% 
#     mutate(Cas9Display = Cas9Display_name) %>% 
#     arrange(n_base)
#   with_rolling_mean <- tmp_mat %>% 
#     mutate(rollingGtoA = zoo::rollapply(GtoA, width = 51, mean, align = "center", na.rm = TRUE, fill = NA)) %>% 
#     mutate(rollingCtoT = zoo::rollapply(CtoT, width = 51, mean, align = "center", na.rm = TRUE, fill = NA)) 
#   
#   my_data[[i]] <- with_rolling_mean
# }
# all_roll_means <- do.call("rbind", my_data)
# 
# # DC573_only <- all_roll_means %>% filter(loci == "HEK3")
# # DC573_only_all_base <- grouped_by_condition %>% filter(loci == "HEK3")
# # ggplot(DC573_only_all_base, aes(n_base, GtoA)) + geom_point() + facet_wrap(~condition)
# 
# 
# ##### Fit Best fit line #####
# data_for_fit_reverse <- all_roll_means %>% filter(n_base >= -1000 & !is.na(helicaseDisplay)) %>% 
#   mutate(logRollingGtoA = log(rollingGtoA)) 
# 
# lm_eqn = function(df){
#   m = lm(logRollingGtoA ~ n_base, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(as.numeric(coef(m)[1]), digits = 2), 
#                         b = format(as.numeric(coef(m)[2]), digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# eq <- ddply(data_for_fit_reverse,.(helicaseDisplay, Cas9Display),lm_eqn)
# 
# 
# ggplot(all_roll_means %>% filter(n_base >= -1000 & !is.na(helicaseDisplay)), aes(n_base, rollingGtoA, color = loci)) + geom_point() + 
#   facet_wrap(~helicaseDisplay+Cas9Display) 
# ggplot(grouped_by_base %>% filter(Cas9Display == "nCas9 D10A"), aes(n_base, meanLogRollingGtoA, color = helicaseDisplay)) +
#   geom_point() + # facet_grid(~Cas9Display) + 
#   scale_color_manual(values = newPalette400_rotate)
# 
# 
# 
# p <- ggplot(data = data_for_fit_reverse, aes(n_base, logRollingGtoA, color = loci)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=FALSE, color="black", formula = y~ x) 
# p1 <- p + geom_text(data=eq,aes(x = -300, y = -6,label=V1), parse = TRUE, inherit.aes=FALSE, size = 5) + 
#   facet_grid(helicaseDisplay~Cas9Display)
# p1full_reverse <- p1 + theme_bw() + scale_color_startrek() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(-9, 0)
# 
# ggsave(file.path(out_dir, "fig2e_range_modelling_reverse.pdf"), p1full_reverse, units = "cm", width = 30, height = 30, dpi = 300)
# 
# 
# g_forward_reverse <- arrangeGrob(p1full_reverse,p1full_forward, nrow = 1)
# ggsave(file.path(out_dir, "fig2e_range_modelling_both.pdf"), g_forward_reverse, units = "cm", width = 60, height = 30, dpi = 300)
# 
# 
# data_for_fit_forward_reverse <- rbind(data_for_fit_reverse, data_for_fit_forward) %>% filter(rollingGtoA > 0)
# 
# data_for_fit_forward_reverse_grouped <- data_for_fit_forward_reverse %>%
#   group_by(n_base, helicaseDisplay, Cas9Display) %>% 
#   dplyr::summarise(n= n(), meanRollingGtoA = mean(rollingGtoA), SD = sd(rollingGtoA), 
#                    meanLogRollingGtoA = mean(logRollingGtoA), logSD = sd(logRollingGtoA))
# 
# 
# p_PcrAM6_only <- ggplot(data = data_for_fit_forward_reverse_grouped %>% filter(helicaseDisplay == "PcrA M6"), aes(n_base, log10(meanRollingGtoA), color = Cas9Display)) +
#   geom_point() +
#   # geom_smooth(method = "lm", se=FALSE, color="black", formula = y~ x) + 
#   # geom_text(data=eq,aes(x = -300, y = -6,label=V1), parse = TRUE, inherit.aes=FALSE, size = 5) + 
#   # facet_wrap(~Cas9Display, nrow = 2) + 
#   theme_bw() + scale_color_startrek() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())+ 
#   xlab("Base Position") + ylab("log(Average Rolling G>A)") + 
#   geom_vline(xintercept=c(1), linetype="dashed")
# 
# ggsave(file.path(out_dir, "fig2e_range_modelling_both_PcrA_M6.pdf"), p_PcrAM6_only, units = "cm", width = 16, height = 8, dpi = 300)
# 
#     # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
# p_PcrAM6_only <- ggplot(data = data_for_fit_forward_reverse %>% filter(helicaseDisplay == "PcrA M6"), aes(n_base, log10(rollingGtoA), color = loci)) +
#   geom_point() +
#   # geom_smooth(method = "lm", se=FALSE, color="black", formula = y~ x) + 
#   # geom_text(data=eq,aes(x = -300, y = -6,label=V1), parse = TRUE, inherit.aes=FALSE, size = 5) + 
#   facet_wrap(~Cas9Display, nrow = 2) + 
#   theme_bw() + scale_color_startrek() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   xlab("Base Position") + ylab("log(Average Rolling G>A)")
#   # geom_smooth(data=dplyr::filter(data_for_fit_forward_reverse,n_base>3), method = 'lm',se = FALSE) # + theme(legend.position = c(0.9, 0.85))  
# ggsave(file.path(out_dir, "fig2e_range_modelling_both_PcrA_M62.pdf"), p_PcrAM6_only, units = "cm", width = 16, height = 16, dpi = 300)


##### OKAY I guess we will try it again this time considering all bases and look at the overall mutation rate. #####
with_mutation_rate_all <- with_plasmid_annotation %>%
  filter((A+C+G+`T`)>5000) %>%
  mutate(A_perc = ifelse(ref_base == 'A', 0, A/base_sum),
         C_perc = ifelse(ref_base == 'C', 0, C/base_sum),
         G_perc = ifelse(ref_base == 'G', 0, G/base_sum),
         T_perc = ifelse(ref_base == 'T', 0, `T`/base_sum),
         GtoA = ifelse(ref_base == "G", A/(A+G), NA),
         CtoT = ifelse(ref_base == "C", `T`/(`T`+C), NA)) %>% 
  mutate(total_perc = A_perc + C_perc + G_perc + T_perc) %>%
  mutate(isControl = ifelse(condition %in% c("RUNX1-Ctrl", "DC203-Control", "DC329-Ctrl", "DC330-Ctrl", "DC331-Ctrl", "DC573-DC573", "DC576-DC573", "DC577-DC573"), 1, 0))

control_rates <- with_mutation_rate_all %>% 
  filter(isControl  ==  1) %>% 
  group_by(chr, n_base, ref_base) %>%
  dplyr::summarise(total_perc_ctrl = mean(total_perc))

with_mut_and_control <- merge(with_mutation_rate_all, control_rates, by = c("chr", "n_base"), all.x = T) %>% 
  as_tibble() %>% 
  arrange(n_base)

top_ctrl <- with_mut_and_control %>% filter(total_perc_ctrl > 0.05)

# Filter out bases where the control rates are more than 5%. 
with_mut_background_filtered <- with_mut_and_control %>% 
  filter(total_perc_ctrl < 0.0025) %>% 
  filter(n_base >= -1500 & n_base <= 1500) 
  # filter(loci %in% c("MAP2K1", "IL6", "RUNX1", "TNF", "HEK3")) 
  # filter(!(loci %in% c("VEGFA", "AKT1", "CD209", "RUNX1")))

grouped_by_condition <- with_mut_background_filtered %>% 
  group_by(condition, n_base, chr, ref_base.y, loci, helicaseDisplay, Cas9Display) %>% 
  dplyr::summarise(editRate = mean(total_perc), sd = sd(total_perc), GtoA = mean(GtoA), CtoT = mean(CtoT)) %>%
  # filter(Cas9Display %in% c("nCas9 D10A","nCas9 H840A")| is.na(Cas9Display)) %>% 
  # filter(n_base >= 7 & n_base < 1000) %>%
  filter(editRate < 0.05) %>% 
  # Scale 100x.
  mutate(editRate = editRate * 100) %>% 
  mutate(GtoA = GtoA *100) %>% 
  mutate(CtoT = CtoT *100) 

ggplot(grouped_by_condition %>% filter(loci == "IL6")) + 
  #geom_point(aes(n_base, editRate), color = "black", alpha = 0.3) + 
  geom_point(aes(n_base, GtoA), color = "#F87171", alpha = 0.9) + 
  geom_point(aes(n_base, CtoT), color = "#60A5FA", alpha = 0.9) + 
  facet_grid(helicaseDisplay~Cas9Display) + # + ylim(0,0.1)
  ggtitle("IL6 (5% max)")

unique_condition <- unique(grouped_by_condition$condition)
all_bases <- data.frame(n_base = seq(0, 1500))
my_data <- list()
for (i in seq_along(unique_condition)){
  print(paste("Current i", i))
  tmp_mat <- grouped_by_condition %>% filter(condition == unique_condition[i]) 
  condition_name = tmp_mat$condition[1]
  loci_name = tmp_mat$loci[1]
  helicaseDisplay_name = tmp_mat$helicaseDisplay[1]
  Cas9Display_name = tmp_mat$Cas9Display[1]
  
  tmp_mat <- merge(all_bases, tmp_mat, by = "n_base", all.x = T) %>% 
    as_tibble() %>%
    filter(n_base < -3 | n_base > 17) %>% 
    mutate(condition = condition_name) %>% 
    mutate(loci = loci_name) %>% 
    mutate(helicaseDisplay = helicaseDisplay_name) %>% 
    mutate(Cas9Display = Cas9Display_name) %>% 
    arrange(n_base)
  with_rolling_mean <- tmp_mat %>% 
    mutate(rollingTotalMut = zoo::rollapply(editRate, width = 51, mean, align = "center", na.rm = TRUE, fill = NA)) 
  my_data[[i]] <- with_rolling_mean
}
all_roll_means <- do.call("rbind", my_data)

DC573_only <- all_roll_means %>% filter(loci == "MAP2K1")
DC573_only_all_base <- grouped_by_condition %>% filter(loci == "MAP2K1")
ggplot(DC573_only_all_base, aes(n_base, editRate)) + geom_point() + facet_wrap(helicaseDisplay~condition, ncol = 3) # + ylim(0,0.1)


##### Fit Best fit line #####
data_for_fit_forward <- all_roll_means %>% filter(n_base <= 800 ) %>% # & !is.na(helicaseDisplay)
  mutate(logRollingTotalMut = log(rollingTotalMut)) 

# If the fit is just on the raw data
# data_for_fit <- grouped_by_condition %>% filter(n_base <= 1000 & !is.na(helicaseDisplay)) %>% 
#   filter(loci != "AKT1") %>% 
#   mutate(logRollingGtoA = log(GtoA)) 


# lm_eqn = function(df){
#   m = lm(logRollingGtoA ~ n_base, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(as.numeric(coef(m)[1]), digits = 2), 
#                         b = format(as.numeric(coef(m)[2]), digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# eq <- ddply(data_for_fit_forward,.(helicaseDisplay, Cas9Display),lm_eqn)


ggplot(all_roll_means %>% filter(n_base <= 800), aes(n_base, log(rollingTotalMut), color = loci)) + 
geom_point() + 
facet_grid(Cas9Display~helicaseDisplay) # + ylim(0,0.01)

# Group by nCas9, helicase, and n_base
grouped_by_base <- data_for_fit_forward %>% 
  dplyr::group_by(helicaseDisplay, Cas9Display, n_base) %>% 
  dplyr::summarise(meanLogRollingTotalMut = log(mean(rollingTotalMut))) %>% 
  dplyr::mutate(log_n_base = log(n_base))

lm_eqn = function(df){
  m = lm(meanLogRollingTotalMut ~ n_base, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(as.numeric(coef(m)[1]), digits = 5), 
                        b = format(as.numeric(coef(m)[2]), digits = 5), 
                        r2 = format(summary(m)$r.squared, digits = 5)))
  as.character(as.expression(eq));                 
}
eq <- ddply(grouped_by_base,.(helicaseDisplay, Cas9Display),lm_eqn)

newPalette400_rotate <- c("#F87171","#60A5FA", "#34D399",  "#FBBF24", "#A3E635", "#E879F9", "#A78BFA")
pD10A <- ggplot(grouped_by_base %>% filter(Cas9Display == "nCas9 D10A"), aes(n_base, meanLogRollingTotalMut, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") # + ylim(-5,0) + theme(legend.position = "none")
pH840A <- ggplot(grouped_by_base %>% filter(Cas9Display == "nCas9 H840A"), aes(n_base, meanLogRollingTotalMut, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") # + ylim(-5,0)+ theme(legend.position = c(0.9, 0.8)) 


g <- arrangeGrob(pD10A, pH840A, nrow = 1)
as_ggplot(g)
# ggsave(file.path(out_dir, "fig2e_range_modelling_helicase_all_mut_modes.pdf"), g, units = "cm", width = 16, height = 6, dpi = 150)


max_first  <- max(data_for_fit_forward$editRate, na.rm = T)   # Specify max of first y axis
max_second <- max(data_for_fit_forward$rollingTotalMut, na.rm = T) # Specify max of second y axis
min_first  <- min(data_for_fit_forward$editRate, na.rm = T)   # Specify min of first y axis
min_second <- min(data_for_fit_forward$rollingTotalMut, na.rm = T) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}


# Plot the rolling average above the points?
ggplot(data_for_fit_forward %>% filter(loci == "RUNX1" & Cas9Display == "nCas9 H840A")) + 
  geom_point(aes(n_base, editRate)) + 
  facet_wrap(~condition, ncol = 1) + 
  geom_line(aes(n_base, inv_scale_function(rollingTotalMut, scale, shift))) + 
  scale_y_continuous(limits = c(min_first, max_first), sec.axis = sec_axis(~scale_function(., scale, shift), name="logRollingTotalMut"))





g2 <- ggplot(grouped_by_base, aes(n_base, meanLogRollingGtoA, color = helicaseDisplay)) +
  geom_point(size = 0.8) + # facet_grid(~Cas9Display) + 
  scale_color_manual(values = newPalette400_rotate) + 
  xlab("Base Position (bp)") + theme(legend.position = "none") + 
  geom_text(data=eq,aes(x = 300, y = -4,label=V1), parse = TRUE, inherit.aes=FALSE, size = 5) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y~ x) + 
  facet_grid(helicaseDisplay~Cas9Display)

ggsave(file.path(out_dir, "fig2e_range_modelling2.pdf"), g2, units = "cm", width = 30, height = 30, dpi = 300)

