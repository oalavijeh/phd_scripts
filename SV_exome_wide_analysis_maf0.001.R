#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly = TRUE)

## This script divides the rare SV calls into types and per desired gene if needed in both the cases and the controls
# outputs a table of p values for further analysis

setwd("/re_gecip/renal/oalavijeh/projects/pkd_sv/analysis/exomewide_nmd")

library(ggplot2)
library(data.table)
library(dplyr)

group="pkd"
bed="exonwide"
nc=275
ncon=25962
svlen=0.05
#pheno=arg[1]

ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}
cc = function(x){length(unique(x[,5]))}

# read input files output of InHouse_MAF script
cases <- read.delim(paste("Rare_StructuralVariants_maf0.001",group,"_filtered_",bed,"_",nc,"x_participants.txt", sep=""),sep="\t")
controls <- read.delim(paste("Rare_StructuralVariants_maf0.001rdcontrol_filtered_",bed,"_",ncon,"x_participants.txt", sep=""),sep="\t")
cases <- na.omit(cases)
controls <- na.omit(controls)

# Remove chrX (only analysing autosomal SVs)
a <- cases[cases$CHROM!="chrX",]
b <- controls[controls$CHROM!="chrX",]

# Rearrange columns and convert Mb to kb
a <- a[,c(8,9,10,4,1,2,3,7,11,5)]
a$SVlengthkb <-an(a$LENGTH.Mb.)*1000
a<- a[,-c(9)]
b <- b[,c(8,9,10,4,1,2,3,7,11,5)]
b$SVlengthkb <- an(b$LENGTH.Mb.)*1000
b<- b[,-c(9)]

# Subset by phenotype
#hpo <- read.table(paste("/re_gecip/renal/mchan/CAKUT/lists/",pheno,".txt", sep=""), header=F)
#names(hpo) <- "Part_ID"
#hpo$Part_ID <- ac(hpo$Part_ID)
#nc=length(unique(hpo$Part_ID))
#a <- inner_join(a,hpo,by="Part_ID")

# Remove common variants seen in cancer and gnomad controls with AF > 0.001. These are SVs merged by survivor within 300bp of start/end points with AC > 24.
# Bedtools intersect with > 70% overlap between all case/control SVs to create a list of SVs that need to be removed, matched by SV type.
cancer <- read.table("pkd_rdcontrol_sv_toremove_cancer_gnomad_0.001.txt", header=F)
cancer <- cancer[cancer$V4!="TRA",]
names(cancer) <- c("CHROM", "STARTPOS", "ENDPOS", "CONSEQUENCE")
filtered_cases <- anti_join(a, cancer, by=c("CHROM", "STARTPOS", "ENDPOS"))
filtered_controls <- anti_join(b, cancer, by=c("CHROM", "STARTPOS", "ENDPOS"))

# Read in list of GENCODE v29 protein coding genes (n=19907)
panel <- read.table("/re_gecip/renal/oalavijeh/projects/pkd_sv/scripts/new_workflowv2_2022/gencode_v29.bed", header=F)
panel <- as.vector(panel$V1)

# Create empty vectors
del_fisher=c()
dup_fisher=c()
ins_fisher=c()
inv_fisher=c()
cnv_fisher=c()
sv_fisher=c()
cnv1=c()
cnv2=c()
cnv3=c()
cnv4=c()
cnv5=c()
case_del=c()
case_dup=c()
control_del=c()
control_dup=c()
med_cnv=c()

# Initialise counter
count=0

for(gene in panel){
  a <- subset(filtered_cases, grepl(paste("\\b",gene,"\\b", sep=""), GeneSymbol))
  b <- subset(filtered_controls, grepl(paste("\\b",gene,"\\b", sep=""), GeneSymbol))

  count=count+1
  print(count)

  #subset into SV type

  del <- subset(a, a$CONSEQUENCE == "DEL" & a$SVlengthkb >= svlen)
  ins <- subset(a, a$CONSEQUENCE == "INS" & a$SVlengthkb >= svlen)
  dup <- subset(a, a$CONSEQUENCE == "DUP" & a$SVlengthkb >= svlen)
  inv <- subset(a, a$CONSEQUENCE == "INV" & a$SVlengthkb >= svlen)
  cnv <- subset(a, a$CONSEQUENCE == "CNV" & a$SVlengthkb >= svlen)
  all <- rbind(del,ins,inv,dup,cnv)

  del_b <- subset(b, b$CONSEQUENCE == "DEL" & b$SVlengthkb >= svlen)
  ins_b <- subset(b, b$CONSEQUENCE == "INS" & b$SVlengthkb >= svlen)
  dup_b <- subset(b, b$CONSEQUENCE == "DUP" & b$SVlengthkb >= svlen)
  inv_b <- subset(b, b$CONSEQUENCE == "INV" & b$SVlengthkb >= svlen)
  cnv_b <- subset(b, b$CONSEQUENCE == "CNV" & b$SVlengthkb >= svlen)
  all_b <- rbind(del_b, ins_b, inv_b, dup_b, cnv_b)

  #create data frame for no individuals with rare (MAF< 0.01) SV >= 50bp and run fishers exact test on each row
  DEL <- c(an(length(unique(del$Part_ID))),nc-an(length(unique(del$Part_ID))),an(length(unique(del_b$Part_ID))), ncon-an(length(unique(del_b$Part_ID))))
  DUP <- c(an(length(unique(dup$Part_ID))),nc-an(length(unique(dup$Part_ID))),an(length(unique(dup_b$Part_ID))), ncon-an(length(unique(dup_b$Part_ID))))
  INS <- c(an(length(unique(ins$Part_ID))),nc-an(length(unique(ins$Part_ID))),an(length(unique(ins_b$Part_ID))), ncon-an(length(unique(ins_b$Part_ID))))
  INV <- c(an(length(unique(inv$Part_ID))),nc-an(length(unique(inv$Part_ID))),an(length(unique(inv_b$Part_ID))), ncon-an(length(unique(inv_b$Part_ID))))
  CNV <- c(an(length(unique(cnv$Part_ID))),nc-an(length(unique(cnv$Part_ID))),an(length(unique(cnv_b$Part_ID))), ncon-an(length(unique(cnv_b$Part_ID))))
  TOTAL <- c(an(length(unique(all$Part_ID))),nc-an(length(unique(all$Part_ID))),an(length(unique(all_b$Part_ID))),ncon-an(length(unique(all_b$Part_ID))))
  d <-t(data.frame(DEL, DUP, INS,INV,CNV, TOTAL))
  colnames(d) <- c("cases_aff", "cases_unaff", "controls_aff", "controls_unaff")
  d <- as.data.frame(d)

  row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
    f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
    return(c(row,
             p_val = f$p.value,
             or = f$estimate[[1]],
             or_ll = f$conf.int[1],
             or_ul = f$conf.int[2]))
  }
  p <- data.frame(t(apply(d, 1, row_fisher)))

  del_fisher <- append(del_fisher, p[1,5])
  dup_fisher <- append(dup_fisher, p[2,5])
  ins_fisher <- append(ins_fisher, p[3,5])
  inv_fisher <- append(inv_fisher, p[4,5])
  cnv_fisher <- append(cnv_fisher, p[5,5])
  sv_fisher <- append(sv_fisher, p[6,5])

  #create bar plot for different SV types vs frequencies in cases and controls

  d$prop <- an(d[,1]/nc*100)
  d$prop_con <- an(d[,3]/ncon*100)

  cohort <- c("case", "control","case", "control","case", "control","case", "control", "case", "control", "case", "control")
  sv <- c("DEL", "DEL", "DUP", "DUP", "INS", "INS", "INV", "INV", "CNV", "CNV", "ALL", "ALL")
  freq <- as.vector(t(cbind(d$prop,d$prop_con)))
  q <- data.frame("COHORT"=cohort, "SV"=sv, "FREQ"=freq)

  # Analyse CNV by size - create data frame
  cnv_size <- c("100-249kb", "250-499kb", "500-999kb", ">=1000kb","ALL")
  cases_aff <- c(cc(cnv[cnv$SVlengthkb >=100 & cnv$SVlengthkb<250,]), cc(cnv[cnv$SVlengthkb >=250 & cnv$SVlengthkb<500,]),cc(cnv[cnv$SVlengthkb >=500 & cnv$SVlengthkb<1000,]),cc(cnv[cnv$SVlengthkb >=1000,]), cc(cnv[cnv$SVlengthkb>=100,]))
  cases_unaff <- c(nc-cc(cnv[cnv$SVlengthkb >=100 & cnv$SVlengthkb<250,]), nc-cc(cnv[cnv$SVlengthkb >=250 & cnv$SVlengthkb<500,]),nc-cc(cnv[cnv$SVlengthkb >=500 & cnv$SVlengthkb<1000,]),nc-cc(cnv[cnv$SVlengthkb >=1000,]),nc-cc(cnv[cnv$SVlengthkb>=100,]))
  controls_aff <-c(cc(cnv_b[cnv_b$SVlengthkb >=100 & cnv_b$SVlengthkb<250,]), cc(cnv_b[cnv_b$SVlengthkb >=250 & cnv_b$SVlengthkb<500,]),cc(cnv_b[cnv_b$SVlengthkb >=500 & cnv_b$SVlengthkb<1000,]),cc(cnv_b[cnv_b$SVlengthkb >=1000,]),cc(cnv_b[cnv_b$SVlengthkb >=100,]))
  controls_unaff <- c(ncon-cc(cnv_b[cnv_b$SVlengthkb >=100 & cnv_b$SVlengthkb<250,]), ncon-cc(cnv_b[cnv_b$SVlengthkb >=250 & cnv_b$SVlengthkb<500,]),ncon-cc(cnv_b[cnv_b$SVlengthkb >=500 & cnv_b$SVlengthkb<1000,]),ncon-cc(cnv_b[cnv_b$SVlengthkb >=1000,]),ncon-cc(cnv_b[cnv_b$SVlengthkb>=100,]))
  f <- data.frame(cnv_size, cases_aff, cases_unaff, controls_aff, controls_unaff)
  rownames(f) <- f[,1]
  f <- f[,-1]
  f <- data.frame(t(apply(f, 1, row_fisher)))

  cnv1<-append(cnv1,f[1,5])
  cnv2<-append(cnv2,f[2,5])
  cnv3<-append(cnv3,f[3,5])
  cnv4<-append(cnv4,f[4,5])
  cnv5<-append(cnv5,f[5,5])

  # Frequency of duplication/deletion for CNV > 100kb

  COHORT <- c("CASE", "CONTROL", "CASE", "CONTROL")
  CNV <- c("DEL", "DEL", "DUP", "DUP")
  FREQ <- c(cc(cnv[cnv$SVlengthkb>=100 & (cnv$GENOTYPE_CN=="CN1" | cnv$GENOTYPE_CN=="CN0"),])/nc*100, cc(cnv_b[cnv_b$SVlengthkb>=100 & (cnv_b$GENOTYPE_CN=="CN1"| cnv_b$GENOTYPE_CN=="CN0"),])/ncon*100, cc(cnv[cnv$SVlengthkb>=100 & (cnv$GENOTYPE_CN=="CN3"| cnv$GENOTYPE_CN=="CN4" | cnv$GENOTYPE_CN=="CN5"),])/nc*100,cc(cnv_b[cnv_b$SVlengthkb>=100 & (cnv_b$GENOTYPE_CN=="CN3"| cnv_b$GENOTYPE_CN=="CN4" | cnv_b$GENOTYPE_CN=="CN5"),])/ncon*100)
  r <- data.frame("COHORT"=COHORT, "CNV"=CNV, "FREQ"=FREQ)

  case_del <- append(case_del, r[1,3])
  case_dup <- append(case_dup, r[3,3])
  control_del <- append(control_del, r[2,3])
  control_dup <- append(control_dup, r[4,3])

  # Wilcoxon test to determine difference in median CNV length
  if (dim(cnv)[1]==0){
    w <- "NA"
  } else if (nrow(cnv[cnv$SVlengthkb>=100,])==0 | nrow(cnv_b[cnv_b$SVlengthkb>=100,])==0) {
    w <- "NA"
  } else {
    cnv$cohort <- "CASE"
    cnv_b$cohort <- "CONTROL"
    cnv_all <- rbind(cnv[cnv$SVlengthkb>=100,c(11,10)], cnv_b[cnv_b$SVlengthkb>=100,c(11,10)])
    w <- wilcox.test(SVlengthkb~cohort, data=cnv_all, exact=FALSE)$p.value
  }
  med_cnv <- append(med_cnv,w)

}

df <- data.frame(panel, del_fisher, dup_fisher, ins_fisher, inv_fisher, cnv_fisher, sv_fisher, cnv1, cnv2, cnv3, cnv4, cnv5, case_del, case_dup, control_del, control_dup, med_cnv)
rownames(df) <- df[,1]
colnames(df) <- c("GENE", "DEL", "DUP", "INS", "INV", "CNV", "ALL_SV", "100-249kb", "250-499kb", "500-999kb", ">=1000kb", "ALL_CNV", "CASE_DEL", "CASE_DUP", "CONTROL_DEL", "CONTROL_DUP", "MEDIAN_CNV_SIZE")
write.table(df, paste(group,"maf<0.001_filtered_exome_wide_SV_summary_stats.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
