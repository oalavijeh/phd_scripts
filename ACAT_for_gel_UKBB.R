library(ACAT)
library(data.table)
library(tidyr)
library(data.table)
library(qqman)
library(sumFREGAT)
gene <- fread("Documents/Projects/stones/reviewer_data/kidney-stone-azphewas-com-450k-phewas-binary.csv")
colnames(gene) <- c("Gene","Phenotype","CollapsingModel","Type","pValue","Category","nSamples","BinNcases","BinQVcases","BinNcontrols","BinQVcontrols","BinCaseFreq","BinCtrlFreq","BinOddsRatio","BinOddsRatioLCI","BinOddsRatioUCI")
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

#metanalysis with Fishers
# load the dplyr package for data manipulation
library(dplyr)

# create a sample data frame
df <- data.frame(
  Gene = c("A", "B", "C", "D", "E"),
  Pvalue = c(0.05, 0.02, 0.01, 0.04, 0.03),
  pValue = c(0.01, 0.03, 0.05, 0.02, 0.04),
  N1 = c(100, 150, 200, 125, 175), # sample size of cohort 1
  N2 = c(150, 200, 250, 175, 225) # sample size of cohort 2
)

# calculate the effect sizes for each P-value
final_table <- final_table %>%
  mutate(
    Z1 = qnorm(1 - Pvalue/2) * sign(-log10(Pvalue)),
    Z2 = qnorm(1 - pValue/2) * sign(-log10(pValue)),
    var1 = 1/N1,
    var2 = 1/N2
  )

# calculate the combined effect size using Fisher's method
final_table <- final_table %>%
  mutate(
    Z = (Z1 * sqrt(var2) + Z2 * sqrt(var1)) / sqrt(var1 + var2),
    P = 2 * (1 - pnorm(abs(Z)))
  )

# print the result
final_table

write.table(final_table, "Documents/Projects/stones/reviewer_data/acat_gene_level_ukbb_gel_fishers.txt", col.names = T, row.names = F,
            quote = F, sep = '\t')
