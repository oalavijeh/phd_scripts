library(tidyverse)
library(ggtext)
library(normentR)
library(data.table)

#gwas <- fread("Documents/Projects/pkd/meta_analysis/hetero_no_overlap_gel_finngenr8_ukbb_jbb_meta_analysis1.tbl")
#gwas <- gwas %>% separate(MarkerName, c("chr","pos","ref","alt"), sep = "_")
#meta$chr <- gsub("X", 23, meta$chr) )
#write.table(gwas, "Documents/Projects/pkd/meta_analysis/hetero_no_overlap_gel_finngenr8_ukbb_jbb_meta_analysis_for_manhattan_fuma.tbl", col.names = T, row.names = F, quote = F, sep = '\t')
gwas <- fread("Documents/Projects/pkd/meta_analysis/hetero_no_overlap_gel_finngenr8_ukbb_jbb_meta_analysis_for_manhattan_fuma.tbl")
#order

chrOrder<-c(paste(1:22,sep=""),"X")
gwas$chr<-factor(gwas$chr, levels=chrOrder)
gwas <- gwas[order(gwas$chr, gwas$pos),]
#cut down data
sig_data <- gwas %>% 
  subset(gwas$`P-value` < 0.05)

notsig_data <- gwas %>% 
  subset(gwas$`P-value` >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)

gwas <- bind_rows(sig_data, notsig_data)
#plot
data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

#axis
axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

#plot
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(p), 
                                  color = as_factor(chr), size = -log10(p))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

print(manhplot)
