#need R/3.6.2 NOT 4

library(tidyr)
library(data.table)
library(qqman)

trait="all_cystic"

setwd("Documents/Projects/pkd/paper/")
x <- fread("../all_cystic-_summary_stats/all_cystic_summary_stats_hg38.txt")
b <- x[,c(1,2,3,14)]

names(b) <- c("CHR", "BP", "SNPID", "P")
b$CHR=gsub("chrX","chr23",b$CHR)
b$CHR=gsub("X",23,b$CHR)

c <- b[!is.na(b$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
c$CHR <- as.numeric(c$CHR)
d <- na.omit(c)
d <- as.data.frame(d)
p=length(unique(d$SNP))

# Plot Manhattan
png("gel_cystic_cohort", type="cairo", units="in", res=300, width=12, height=6)
par(omi=c(0,0.5,0,0))
manhattan(d, ylim=c(0,10), cex.axis=0.6, col = c("darkorange", "deepskyblue3"), snp = "SNPID",genomewideline =-log10(0.05/p), suggestiveline = FALSE, chrlabs = c(1:22, "X"))
dev.off()


###ggplot
library(ggplot2)
library(ggmanh)
library(tidyverse)
library(data.table)
library(ggtext)
library(ggbreak)
library(GWASTools)

#prep
data_cum <- d %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

d <- d %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add )
#set plotting axes
axis_set <- d%>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- d%>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)

#0.05/length(summary stats)
sig <- 0.05/p

#plotting
manhplot <- ggplot(d, aes(x = BP, y = -log10(d$P), 
                             color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(P-value)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot 
t <- manhplot + scale_y_cut(breaks=c(0,10)) 
t
y <- t + scale_y_cut(breaks=c(149,150)) 
y
o<- y + scale_y_cut(breaks=c(300,310), space = 0.2)
o

#qq
#Genomic inflation factor function
lambda<-function(pvalues){
  chisq <- qchisq(1-pvalues,1)
  lambda=median(chisq,na.rm=TRUE)/qchisq(0.5,1)
  return(lambda)
}
png("metal_metanalysis_qqplot.png", type = "cairo", units="in", res=300, width=6, height=6)
par(omi=c(0,0.5,0,0))
qqPlot(gwas$`P-value`,main=paste0("lambda=",round(lambda(gwas$`P-value`),4)), cex=0.6)
dev.off()

qq(gwas$`P-value`)



