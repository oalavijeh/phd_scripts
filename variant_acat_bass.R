library(ACAT)
library(data.table)
library(tidyr)
library(data.table)
library(qqman)
library(sumFREGAT)
library(MetaSTAAR)
library(GenomicRanges)
library(rtracklayer)
setwd("Documents/Projects/stones/reviewer_data")

bass1 <- fread("bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_missense_nvar3864295_low0_high1e-05_group1.txt")
bass2 <- fread("bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_missense_nvar934162_low1e-05_high0.0001_group2.txt")
bass3 <- fread("bhr/~/dec_bhr_ms_variant_ss_400k_final_thin_withnullburden_missense_nvar130359_low0.0001_high0.001_group3.txt")

bass <- rbind(bass1,bass2,bass3)
bass_acat <- bass[,c("markerID","Pvalue","gene","BETA","AF","annotation")]
colnames(bass_acat) <- c("ID","P","GENE","BETA","AF","annotation")
#sort columns
bass_acat <- separate(bass_acat, "ID", into = c("CHROM","POS"), sep = ":", remove = FALSE)
bass_acat <- separate(bass_acat, "POS", into = c("REF","EA"), sep = "/", remove = FALSE)
bass_acat$REF <- gsub("^.*_","",bass_acat$REF)
bass_acat$POS <- gsub("[^0-9]","",bass_acat$POS)
#order
chrOrder<-c(paste("chr",1:22,sep=""),"chrX", "chrY")
bass_acat$CHROM<-factor(bass_acat$CHROM, levels=chrOrder)
bass_acat$CHROM
bass_acat <- bass_acat[order(bass_acat$CHROM, bass_acat$POS),]

#read in sumfregate bild 38 file 
ref <- fread("build38_european_1000G_sumfregat_reference.txt")

#prepfile
prep.score.files(bass_acat, output.file.prefix = "test", reference = ref)


#####test
# split the data frame by gene
# create a list of p-values per gene
test_acat = bass_acat
test_acat$P <- ifelse(bass_acat$P==1, 0.99,bass_acat$P)
gene_pvals_list <- split(test_acat$P, test_acat$GENE)

# define a function that applies ACAT to a list of p-values
acat_gene <- function(p_values) {
  return(ACAT::ACAT(Pvals = p_values))
}

# apply the function to each element of the list
gene_acat_list <- lapply(gene_pvals_list, acat_gene)

# convert the results to a data frame
gene_acat_df <- data.frame(gene = names(gene_acat_list),
                           p_value = unlist(gene_acat_list))
