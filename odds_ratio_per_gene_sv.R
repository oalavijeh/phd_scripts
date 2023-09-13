#script to work out OR for SV hits 
setwd("~/re_gecip/renal/oalavijeh/projects/pkd_sv/analysis/exomewide/maf0.001/combined_filtered/")
library(ggplot2)
library(data.table)
library(dplyr)

ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}
cc = function(x){length(unique(x[,5]))}

# read input files output of InHouse_MAF script
cases <- fread("Rare_StructuralVariants_maf0.001pkd_filtered_exonwide_1210x_participants.txt")
controls <- fread("Rare_StructuralVariants_maf0.001combined_exonwide_27173x_participants.txt")
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

# Remove common variants seen in cancer and gnomad controls with AF > 0.001. These are SVs merged by survivor within 300bp of start/end points with AC > 24.
# Bedtools intersect with > 70% overlap between all case/control SVs to create a list of SVs that need to be removed, matched by SV type.
cancer <- read.table("../../../exomewide_nmd/pkd_rdcontrol_sv_toremove_cancer_gnomad_0.001.txt", header=F)
cancer <- cancer[cancer$V4!="TRA",]
names(cancer) <- c("CHROM", "STARTPOS", "ENDPOS", "CONSEQUENCE")
filtered_cases <- anti_join(a, cancer, by=c("CHROM", "STARTPOS", "ENDPOS"))
filtered_controls <- anti_join(b, cancer, by=c("CHROM", "STARTPOS", "ENDPOS"))

write.table(filtered_cases, "filtered_cases_all_cystic_for_or.txt", col.names = T, row.names = F, quote = F, sep = '\t')
write.table(filtered_controls, "filtered_controls_all_cystic_for_or.txt", col.names = T, row.names = F, quote = F, sep = '\t')
#all_Sv
gene="HNF1B"
count_case <- filtered_cases[grepl(gene,filtered_cases$GeneSymbol),]
count_control <- filtered_controls[grepl(gene,filtered_controls$GeneSymbol),]
#dels
gene="PKD1"
count_case <- filtered_cases[grepl(gene,filtered_cases$GeneSymbol),]
count_case <- subset(count_case, count_case$CONSEQUENCE=="DEL")
count_control <- filtered_controls[grepl(gene,filtered_controls$GeneSymbol),]
count_control <- subset(count_control, count_control$CONSEQUENCE=="DEL")
#CNV
gene="HNF1B"
type="CNV"
count_case <- filtered_cases[grepl(gene,filtered_cases$GeneSymbol),]
count_case <- count_case[grepl(type,count_case$CONSEQUENCE),]
count_control <- filtered_controls[grepl(gene,filtered_controls$GeneSymbol),]
count_control <- count_control[grepl(type,count_control$CONSEQUENCE),]

library(epitools)
tag <- c("case","control")
outcome <- c("gene_variant","no_gene_variant")
case_n_total <- 1209
control_n_total <- 26096
case_event <- length(unique(count_case$Part_ID)) 
control_event <- length(unique(count_control$Part_ID))
case_nonevent <- (case_n_total-case_event)
control_nonevent <- (control_n_total-control_event)
data<- matrix(c(case_event,case_nonevent,control_event,control_nonevent),nrow=2,ncol=2,byrow=TRUE)
dimnames(data) <- list('tag'=tag, 'Outcome'=outcome)
data
oddsratio(data)

#if small numbers can do fishers
fisher.test(data)

