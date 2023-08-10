library(tidyverse)
library("AnnotationDbi")
library("org.Hs.eg.db")
#setwd
setwd("~/Documents/Projects/stones/reviewer_data/bhr/")
#read in baseline/annotation files
baseline_model <- read.table("~/Documents/Projects/stones/reviewer_data/bhr/ms_baseline_oe5.txt")
annotatoin_1 <- fread("annotation1.txt")
#liability
#calculate sample prevalence
bp_sample_prevalence = 3147/(3147 + 255496) 

#this population prevalence estimate is from Ferrari et al, 2016 Bipolar Disorders
population_prevalence_bp = 0.1 


obs2lia_factor <- function(K, P){
  X <- qnorm(K,lower.tail=FALSE)
  z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
  factor <- (K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
  return(factor)
}

bp_scalingfactor = obs2lia_factor(population_prevalence_bp,bp_sample_prevalence)
#read in summary stats from AZ
az <- fread("../kidney-stone-variant-level.csv")
colnames(az) <- c("Variant","Variant type","Phenotype","Category","Model","Consequence type","Gene",
                     "Transcript","cDNA change","Amino acid change","Exon rank","No. cases","No. AA cases",
                     "No. AB cases","No. BB cases","Case MAF","% AB or BB cases","% BB cases","No. controls",
                     "No. AA controls","No. AB controls","No. BB controls","Control MAF","% AB or BB controls",
                     "% BB controls","p-value","Odds ratio","Odds ratio LCI","Odds ratio UCI")
#transform consequence type

to_repalce<- c("5_prime_UTR_premature_start_codon_gain_variant", 
                "disruptive_inframe_deletion",
                 "frameshift_variant",
                 "splice_acceptor_variant",
                 "splice_donor_variant",
                 "start_lost",
                 "stop_gained",
                "stop_lost")
new_phrase <- "ptv"

az$`Consequence type`<- ifelse(az$`Consequence type` %in% to_repalce, new_phrase, az$`Consequence type`)
table(az$`Consequence type`)
wrangle_sumstats <- function(table,n_cases,n_controls, var_filter) {
  
  #Filter to variants of interest
  table = table[var_filter,] 
  table$ac_case = table$`No. AB cases`+table$`No. BB cases`
  table$ac_ctrl = table$`No. AB controls`+table$`No. BB controls`
    #Compute sample prevalence, will be used to compute per-sd beta
  prevalence = n_cases/(n_cases+n_controls) 
  
  #Compute variant MAF in cases, will be used to compute per-sd beta
  table$AF_case = table$`Case MAF` 
  #Compute variant MAFoverall, will be used to compute per-sd beta
  table$AF = (table$ac_case + table$ac_ctrl)/(2*(n_cases + n_controls)) 
  
  #calculate per-sd betas
  table$beta = (2 * (table$AF_case - table$AF) * prevalence)/sqrt(2 * table$AF * (1 - table$AF) * prevalence * (1 - prevalence))  
  
  #calculate variant variances
  table$twopq = 2*table$AF * (1 - table$AF) 
  
  #convert betas from per-sd (i.e. sqrt(variance explained)) to per-allele (i.e. in units of phenotype) betas.
  #per-allele are the usual betas reported by an exome wide association study.
  table$beta_perallele = table$beta/sqrt(table$twopq)
  #get rid of sex chr 
  table <- table[grepl("^\\d", table$Variant),]
  
  #aggregate into gene-level table.
  #position doesn't need to be super precise, as it is only used to order genes for jackknife
  #N = sum of case and control counts
  sumstats = data.frame(gene = gsub("'","",table$Gene),
                        AF = table$AF,
                        beta = table$beta_perallele,
                        gene_position = parse_number(sapply(strsplit(table$Variant, split = "-"), function(x) x[[2]])),
                        chromosome =  parse_number(sapply(strsplit(table$Variant, split = "-"), function(x) x[[1]])),
                        N = n_cases + n_controls,
                        phenotype_key = "USD")
  
  #we have found that in these smaller sample analyses, there are some genes with
  #large burden scores that are clearly outliers
  #we remove one such gene, TTN
  sumstats$ensid = mapIds(org.Hs.eg.db,
                       keys=sumstats$gene, 
                       column="ENSEMBL",
                       keytype="SYMBOL",
                       multiVals="first")
  
  sumstats$gene <- sumstats$ensid
  sumstats <- subset(sumstats, select=-c(ensid))
  sumstats <- na.omit(sumstats)
  

  }

az_ptv <- wrangle_sumstats(az, 3147,255496, az$`Consequence type`=='ptv')
az_missense <- wrangle_sumstats(az, 3147,255496, az$`Consequence type`=='missense_variant')
#split table into MAF bins 
az_ptv1<- subset(az_ptv, az_ptv$AF>0 & az_ptv$AF<1e-05)
az_ptv2<- subset(az_ptv, az_ptv$AF>1e-05 & az_ptv$AF<1e-04)
az_ptv3<- subset(az_ptv, az_ptv$AF>1e-04 & az_ptv$AF<1e-03)
az_ptv4<- subset(az_ptv, az_ptv$AF>1e-03 & az_ptv$AF<1e-02)

az_missense1<- subset(az_missense, az_missense$AF>0 & az_missense$AF<1e-05)
az_missense2<- subset(az_missense, az_missense$AF>1e-05 & az_missense$AF<1e-04)
az_missense3<- subset(az_missense, az_missense$AF>1e-04 & az_missense$AF<1e-03)
az_missense4<- subset(az_missense, az_missense$AF>1e-03 & az_missense$AF<1e-02)

#put into lists
az_ptv_list <- list(az_ptv2,az_ptv3,az_ptv4)
az_missense_list <- list(az_missense2,az_missense3,az_missense4)

ptv_list_unique <- lapply(az_ptv_list, function(df) df[!duplicated(df), ])
missense_list_unique <- lapply(az_missense_list, function(df) df[!duplicated(df), ])

for (i in seq_along(ptv_list_unique)) {
  assign(paste0("ptv", i), ptv_list_unique[[i]])
}

for (i in seq_along(missense_list_unique)) {
  assign(paste0("missense", i), missense_list_unique[[i]])
}

#run bhr in single cohorts  first
usd_ptv_bhr1 <- BHR(trait1_sumstats = ptv1, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2

usd_ptv_bhr2 <- BHR(trait1_sumstats = ptv2, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2

usd_ptv_bhr3 <- BHR(trait1_sumstats = az_ptv4, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2

usd_missense_bhr1 <- BHR(trait1_sumstats = az_missense2, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2


usd_all_ptv <- BHR(mode = "aggregate", 
                   ss_list = list(az_ptv2,az_ptv3,az_ptv4),
                   trait_list = list("USD"),
                   annotations = list(baseline_model),
                   genomewide_correction = TRUE)

usd_all_ptv
#convert observed scale to liability scale h2

#run BHR on a single fraction 
usd_ptv_bhr1 <- BHR(trait1_sumstats = az_ptv, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2

#run BHR aggregate 
usd_all_ptv <- BHR(mode = "aggregate", 
                   ss_list = table_list,
                   trait_list = list("USD"),
                   annotations = list(baseline_model),
                   genomewide_correction = TRUE)

usd_all_missense <- BHR(mode = "aggregate", 
                   ss_list = table_list,
                   trait_list = list("USD"),
                   annotations = list(baseline_model),
                   genomewide_correction = TRUE)

#liability adjusted h2 across all PTVs 
usd_all_ptv$aggregated_mixed_model_h2* bp_scalingfactor
#associated SE adjusted 
usd_all_ptv$aggregated_mixed_model_h2se^bp_scalingfactor


print(paste0("USD total PTV Burden Heritability: ",
             round(usd_all_ptv$aggregated_mixed_model_h2* bp_scalingfactor, 4),
             " ",
             "(",
             round(usd_all_ptv$aggregated_mixed_model_h2se * bp_scalingfactor, 4),
             ")"))

#with SLC34A3(ENSG00000198569) as a fixed gene effect 
usd_all_ptv <- BHR(mode = "aggregate", 
                   ss_list = list(output_table1,output_table2,output_table3),
                   trait_list = list("USD"),
                   annotations = list(baseline_model),
                   genomewide_correction = TRUE, 
                   fixed_genes = c("ENSG00000198569"))

#bhr per MAF pLOF bin to look for major gene 

usd_ptv_bhr1 <- BHR(trait1_sumstats = output_table1, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2


usd_ptv_bhr2 <- BHR(trait1_sumstats = output_table2, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2

usd_ptv_bhr3 <- BHR(trait1_sumstats = output_table3, 
                    annotations = list(baseline_model), #baseline model including constraint annotations
                    num_blocks = 100, #number of blocks for jackknife
                    mode = "univariate") #run in univariate mode to compute burden h2



