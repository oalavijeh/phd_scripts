library(ACAT)
library(data.table)
library(tidyr)
library(data.table)
library(qqman)
library(sumFREGAT)
library(MetaSTAAR)
setwd("Documents/Projects/stones/reviewer_data")
ukbbv <- fread("kidney-stone-variant-level.csv")
colnames(ukbbv) <- c("Variant","Variant type","Phenotype","Category","Model","Consequence type","Gene",
                    "Transcript","cDNA change","Amino acid change","Exon rank","No. cases","No. AA cases",
                    "No. AB cases","No. BB cases","Case MAF","% AB or BB cases","% BB cases","No. controls",
                    "No. AA controls","No. AB controls","No. BB controls","Control MAF","% AB or BB controls",
                    "% BB controls","p-value","Odds ratio","Odds ratio LCI","Odds ratio UCI")


gene_ptvraredmg <- subset(gene, gene$CollapsingModel=="flexnonsynmtr")
gel <- fread("Documents/Projects/stones/saige_missense_PHENO_hes_stones_maf0.001_missense_cadd_output.txt")
gel <- separate(gel, col="Gene", into=c("Gene", "ENS"), sep="\\_")

ukbb_gel <- merge(gel, gene_ptvraredmg, by.x = "Gene", by.y = "Gene")
ukbb_gel <- ukbb_gel[,c(1,3,26)]
#need to sort out P=1 and then try an apply function 
ACAT(Pvals = ukbb_gel$Pvalue)
ukbb_gel$Pvalue <- ifelse(ukbb_gel$Pvalue == 1.0000000, 0.99, ukbb_gel$Pvalue)

#function to pass
acat_row <- function(row) {
  ACATO(as.numeric(row[2:3]))
}

# Apply the function to each row using apply()
result <- apply(ukbb_gel, 1, acat_row)
final_table <- cbind(ukbb_gel, result)

write.table(final_table, "Documents/Projects/stones/reviewer_data/acat_gene_level_ukbb_gel.txt", col.names = T, row.names = F,
            quote = F, sep = '\t')

