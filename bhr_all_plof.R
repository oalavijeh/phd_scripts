library(tidyverse)
library("AnnotationDbi")
library("org.Hs.eg.db")

setwd("~/Documents/Projects/stones/reviewer_data/bhr/")

baseline_model <- read.table("~/Documents/Projects/stones/reviewer_data/bhr/ms_baseline_oe5.txt")
annotatoin_1 <- fread("annotation1.txt")
plof_variants1 <- bigreadr::fread2("~/Documents/Projects/stones/reviewer_data/bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_pLoF_nvar427542_low0_high1e-05_group1.txt")
plof_variants2 <- bigreadr::fread2("~/Documents/Projects/stones/reviewer_data/bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_pLoF_nvar58328_low1e-05_high0.0001_group2.txt")
plof_variants3 <- bigreadr::fread2("~/Documents/Projects/stones/reviewer_data/bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_pLoF_nvar6215_low0.0001_high0.001_group3.txt")
#save tables
table_list <- list(plof_variants1,plof_variants2, plof_variants3)

#write function
wrangle_genebass <- function(table){
  table <- subset(table, table$phenocode=="132036")
  table = table[,c(4,33,43,1,7,8,5,27)]
  table <- table[grepl("^chr\\d", table$locus),]
  #make table
  table <- data.frame(gene = table$gene,
                         AF = table$AF,
                         beta = table$beta_per_allele,
                         gene_position = parse_number(sapply(strsplit(table$locus, split = ":"), function(x) x[[2]])),
                         chromosome =  parse_number(sapply(strsplit(table$locus, split = ":"), function(x) x[[1]])),
                         N = table$n_cases+table$n_controls,
                         phenotype_key = "USD")
  
  #get ENSEMBL IDS
  table$ensid = mapIds(org.Hs.eg.db,
                          keys=table$gene, 
                          column="ENSEMBL",
                          keytype="SYMBOL",
                          multiVals="first")
  
  table$gene <- table$ensid
  table <- subset(table, select=-c(ensid))
  table <- na.omit(table)
} 


# Create an empty list to store the output tables
output_tables <- list()

# Loop over the input tables and apply the wrangle_genebass function to each
for (i in 1:length(table_list)) {
  output_tables[[i]] <- wrangle_genebass(table_list[[i]])
}

# Assign each output table to a new variable name
output_table1 <- output_tables[[1]]
output_table2 <- output_tables[[2]]
output_table3 <- output_tables[[3]]

#run BHR as an aggregate across all MAF bins

usd_all_ptv <- BHR(mode = "aggregate", 
    ss_list = list(output_table1,output_table2,output_table3),
    trait_list = list("USD"),
    annotations = list(baseline_model),
    genomewide_correction = TRUE)

usd_all_ptv
#convert observed scale to liability scale h2

#calculate sample prevalence
bp_sample_prevalence = 5879/(5879 + 388962) 

#this population prevalence estimate is from Ferrari et al, 2016 Bipolar Disorders
population_prevalence_bp = 0.1 


obs2lia_factor <- function(K, P){
  X <- qnorm(K,lower.tail=FALSE)
  z <- (1/sqrt(2*pi))*(exp(-(X**2)/2))
  factor <- (K*(1-K)*K*(1-K))/(P*(1-P)*(z**2))
  return(factor)
}

bp_scalingfactor = obs2lia_factor(population_prevalence_bp,bp_sample_prevalence)
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



