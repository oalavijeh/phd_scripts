setwd("pkd/meta_analysis/second_attempt_wth_adam_script")
library(data.table)
library(qqman)
d=fread("../finngenr8_cystic_b38_ss.txt",data.table=F)

i="FinnGen"

#d = d[,c("chr","pos","chr_pos","pval")]
#names(d) = c("CHR","BP","SNP","P")
d=d[which(d$chrom %in%c(1:22)),]
d$pos=as.numeric(as.character(d$pos))
d$BP=as.numeric(as.character(d$BP))
d$pval=as.numeric(as.character(d$pval))
names(d)= c("CHR","BP","REF","ALT","RSID","NEAREST_GENES","P","MLOGP","BETA","SEBETA","AF_ALT","AF_ALT_CASES"
            ,"AF_ALT_CONTROLS","SNP")

observed <- sort(d$`P-value`)
lobs <- -(log10(observed))
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
obs.chisq = qchisq(1-d$`P-value`,1)
lambda = round(median(obs.chisq)/qchisq(0.5,1),3)

png(paste(i,"_QQ.png",sep=""),width=6,height=6, units="in",res=300)
par(omi=c(0,0.5,0,0))
plot(lexp,lobs,main=paste("Lambda",lambda,sep="="),xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",cex.axis=0.8)
abline(a=0,b=1)
dev.off()

png(paste(i,"_Manhattan.png",sep=""),units="in",width=12,height=6,res=300)
par(omi=c(0,0.5,0,0))
manhattan(d,suggestiveline = F,cex=0.6,cex.axis=0.8,genomewideline=-log10(5e-8))
dev.off()

#find common markers 
finngen <- fread("finngenr8_cystic_b38_ss_for_matching.txt")
finngen <- finngen$SNP_POS

gel <- fread("gel_b38_cystic_ss_for_matching.tsv")
gel <- gel$SNP_POS

jbb_ukbb <- fread("jbb_ukbb_b38_cystic_ss_for_matching.tsv")
jbb_ukbb <- jbb_ukbb$SNP_POS 

#get common snps 
common_snps <- Reduce(intersect, list(finngen, gel, jbb_ukbb))
#write table
write.table(as.data.frame(common_snps), "common_snps.txt", col.names = F, row.names = F, quote = F)
  
#make each files just the common snps between them 
finngen <- subset(finngen, finngen$SNP_POS %in% common_snps)
gel <- subset(gel, gel$SNP_POS %in% common_snps)
jbb_ukbb <- subset(jbb_ukbb, jbb_ukbb$SNP_POS %in% common_snps)

#write these out and do metanalysis
write.table(finngen, "finngenr8_cystic_b38_commons_snps.txt", col.names = T, row.names = F, quote = F, sep = '\t')
write.table(gel, "gel_cystic_b38_commons_snps.txt", col.names = T, row.names = F, quote = F, sep = '\t')
write.table(jbb_ukbb, "jbb_ukbb_cystic_b38_commons_snps.txt", col.names = T, row.names = F, quote = F, sep = '\t')
#read in meta file 
meta <- fread("stderr_common_snps_metanalysis1.tbl")
#qc
bad = subset(meta, meta$MaxFreq > meta$MinFreq + 0.075 | meta$HetPVal<1e-5)
#remove bad markers
gwas <- subset(meta, !meta$MarkerName %in% bad$MarkerName)
#plot markers
gwas <- gwas %>% separate(MarkerName, c("chr","pos","ref","alt"), sep = "_")
#gwas$chr <- gsub("X", 23, meta$chr) 
#order
chrOrder<-c(paste(1:22,sep=""),"X")
gwas$chr<-factor(gwas$chr, levels=chrOrder)
gwas <- gwas[order(gwas$chr, gwas$pos),]
#save
write.table(gwas, "stderr_common_snps_metanalysis1_good_markers_sorted_b38_cystic.tbl", col.names = T, row.names = F, quote = F, sep = '\t')
gwas$pos <- as.numeric(gwas$pos)
#cut down data
sig_data <- gwas %>% 
  subset(gwas$`P-value` < 0.05)

notsig_data <- gwas %>% 
  subset(gwas$`P-value` >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)

gwas <- bind_rows(sig_data, notsig_data)
#plot
data_cum <- gwas %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas <- gwas %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = pos + bp_add)

#axis
axis_set <- gwas %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas %>% 
  filter(p==min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

#plot
manhplot <- ggplot(gwas, aes(x = bp_cum, y = -log10(p), 
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

#plot using qqman 
meta <- fread("stderr_common_snps_metanalysis1_good_markers_sorted_b38_cystic.tbl")
png("metanalysis_manhattan_qqmanv2.png", type="cairo", units="in", res=300, width=12, height=6) 
par(omi=c(0,0.5,0,0))
manhattan(meta, chr="chr", bp="pos", p= "P-value", snp="snp", suggestiveline = FALSE, 
          genomewideline = -log10(0.00000005), ylim=c(0,8),cex.axis=0.6, cex=0.6, col=c("darkorange", "deepskyblue3"))
dev.off()


ukbb <- fread("jbb_ukbb_b38_cystic_ss_for_matching.tsv")
ukbb$hm_chrom <- ifelse(ukbb$hm_chrom=="X",23,ukbb$hm_chrom)
ukbb <- na.omit(ukbb)
png("ukbb_jbb_manhattan_qqman.png", type="cairo", units="in", res=300, width=12, height=6) 
par(omi=c(0,0.5,0,0))
manhattan(ukbb, chr="chromosome", bp="base_pair_location", p= "p_value", snp="SNP", suggestiveline = FALSE, 
          genomewideline = -log10(0.00000005), ylim=c(0,8),cex.axis=0.6, cex=0.6, col=c("darkorange", "deepskyblue3"))
dev.off()

png("ukbb_jbb_manhattan_qqman.png", type="cairo", units="in", res=300, width=12, height=6) 
par(omi=c(0,0.5,0,0))
qq(ukbb$p_value)
dev.off()

#inflation
chisq <- qchisq(1-ukbb$p_value,1)
median(chisq)/qchisq(0.5,1)


#gel
gel <- fread("gel_b38_cystic_ss_for_matching.tsv")
gel$CHR <- ifelse(gel$CHR=="X",23,gel$CHR)
gel$CHR <- as.numeric(gel$CHR)
png("gel_manhattan_qqman.png", type="cairo", units="in", res=300, width=12, height=6) 
par(omi=c(0,0.5,0,0))
manhattan(gel, chr="CHR", bp="POS", p= "p.value", snp="SNP", suggestiveline = FALSE, 
          genomewideline = -log10(0.00000005), ylim=c(0,8),cex.axis=0.6, cex=0.6, col=c("darkorange", "deepskyblue3")
          ,chrlabs = c(1:22,"X"))
dev.off()

png("ukbb_jbb_manhattan_qqman.png", type="cairo", units="in", res=300, width=12, height=6) 
par(omi=c(0,0.5,0,0))
qq(ukbb$p_value)
dev.off()

#inflation
chisq <- qchisq(1-ukbb$p_value,1)
median(chisq)/qchisq(0.5,1)
