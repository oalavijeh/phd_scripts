#script to take the final_hes table with age of ESRF and combine with a table of worst gene consequence from the tiereing data. 
#first take table of all tiering data for proband cystic kidney disease patients as per pkd_tiering.R <- f
#take all the age of ESRF data <-g as per all_CKD_stage.R

#long method - PKD1 truncating, PKD2 truncating, PKD1 missense, PKD2 missense, PKD1/2 other, other gene (tier1), NMD 
# -	PKD1 Truncating: stop_lost/stop_gained/splice_region_variant/splice_donor_variant/ splice_acceptor_variant /frameshift_variant
# -	PKD2 Truncating: stop_lost/stop_gained/splice_region_variant/splice_donor_variant/ splice_acceptor_variant /frameshift_variant
# -	PKD1+2Missense: self-explanatory 
# -	PKD1/2 other:  all the others below if in PKD1 or 2
# -	Other genes: self-explanatory 
# -	NMD: no TIER1 variants or no PKD1/2 variants 

#read in stuff
library(data.table)
f <- fread("phenotypes/outputs/pkd_tiering_output_f_v16.txt") 
#g<- fread("phenotypes/outputs/all_CKD_stages_output_g.txt") 
final_hes_combined <- fread("phenotypes/outputs/all_races_CKD_stages_from_hes_v16.txt")

# -	Truncating/protein truncating PKD1 and 2: stop_lost/stop_gained/splice_region_variant/splice_donor_variant/ splice_acceptor_variant /frameshift_variant
pkd1_truncating <- subset(f, tier == 'TIER1' & genomic_feature_hgnc == 'PKD1' & (consequence_type == 'stop_gained' | 
                                                                consequence_type == 'stop_lost'|
                                                                consequence_type == 'splice_region_variant'| 
                                                                consequence_type == 'incomplete_terminal_codon_variant'|
                                                                consequence_type == 'start_lost'|
                                                                consequence_type == 'feature_truncation'| 
                                                                consequence_type == 'splice_donor_variant' |
                                                                consequence_type == 'splice_acceptor_variant' |
                                                                consequence_type == 'frameshift_variant'))
pkd1_truncating <- pkd1_truncating[!duplicated(pkd1_truncating$participant_id),]
pkd1_truncating$type <- rep('pkd1_truncating', nrow(pkd1_truncating))

pkd2_truncating <- subset(f, tier == 'TIER1' & genomic_feature_hgnc == 'PKD2' & (consequence_type == 'stop_gained' | 
                                                                                   consequence_type == 'stop_lost'|
                                                                                   consequence_type == 'splice_region_variant'| 
                                                                                   consequence_type == 'incomplete_terminal_codon_variant'|
                                                                                   consequence_type == 'start_lost'|
                                                                                   consequence_type == 'feature_truncation'| 
                                                                                   consequence_type == 'splice_donor_variant' |
                                                                                   consequence_type == 'splice_acceptor_variant' |
                                                                                   consequence_type == 'frameshift_variant'))
pkd2_truncating <- pkd2_truncating[!duplicated(pkd2_truncating$participant_id),]
pkd2_truncating$type <- rep('pkd2_truncating', nrow(pkd2_truncating))

# -	Non-protein truncating: Missense: missense_variant 
pkd1_non_truncating <- subset(f, tier == 'TIER1' & genomic_feature_hgnc == 'PKD1' & (consequence_type == 'missense_variant' | 
                                                                                       consequence_type == 'coding_sequence_variant' |
                                                                                       consequence_type == 'intron_variant' |
                                                                                       consequence_type == 'non_coding_transcript_exon_variant' |
                                                                                       consequence_type == 'upstream_variant' |
                                                                                       consequence_type == '2KB_upstream_variant' |
                                                                                       consequence_type == 'downstream_gene_variant' |
                                                                                       consequence_type == 'incomplete_terminal_codon_variant' |
                                                                                       consequence_type == 'non_coding_transcript_variant' |
                                                                                       consequence_type == 'stop_retained_variant' |
                                                                                       consequence_type == '2KB_downstream_variant' |
                                                                                       consequence_type == '3_prime_UTR_variant' |
                                                                                       consequence_type == 'downstream_variant' |
                                                                                       consequence_type == 'inframe_deletion' |
                                                                                       consequence_type == 'start_retained_variant' |
                                                                                       consequence_type == 'synonymous_variant' |
                                                                                       consequence_type == '5_prime_UTR_variant' |
                                                                                       consequence_type == 'inframe_insertion' |
                                                                                       consequence_type == 'NMD_transcript_variant' |
                                                                                       consequence_type == 'inframe_insertion' |
                                                                                       consequence_type == 'upstream_gene_variant'))

pkd1_non_truncating <- pkd1_non_truncating[!duplicated(pkd1_non_truncating$participant_id),]
pkd1_non_truncating$type <- rep('pkd1_non_truncating', nrow(pkd1_non_truncating))

pkd2_non_truncating <- subset(f, tier == 'TIER1' & genomic_feature_hgnc == 'PKD2' & (consequence_type == 'missense_variant' | 
                                                                                       consequence_type == 'coding_sequence_variant' |
                                                                                       consequence_type == 'intron_variant' |
                                                                                       consequence_type == 'non_coding_transcript_exon_variant' |
                                                                                       consequence_type == 'upstream_variant' |
                                                                                       consequence_type == '2KB_upstream_variant' |
                                                                                       consequence_type == 'downstream_gene_variant' |
                                                                                       consequence_type == 'incomplete_terminal_codon_variant' |
                                                                                       consequence_type == 'non_coding_transcript_variant' |
                                                                                       consequence_type == 'stop_retained_variant' |
                                                                                       consequence_type == '2KB_downstream_variant' |
                                                                                       consequence_type == '3_prime_UTR_variant' |
                                                                                       consequence_type == 'downstream_variant' |
                                                                                       consequence_type == 'inframe_deletion' |
                                                                                       consequence_type == 'start_retained_variant' |
                                                                                       consequence_type == 'synonymous_variant' |
                                                                                       consequence_type == '5_prime_UTR_variant' |
                                                                                       consequence_type == 'inframe_insertion' |
                                                                                       consequence_type == 'NMD_transcript_variant' |
                                                                                       consequence_type == 'inframe_insertion' |
                                                                                       consequence_type == 'upstream_gene_variant'))

pkd2_non_truncating <- pkd2_non_truncating[!duplicated(pkd2_non_truncating$participant_id),]
pkd2_non_truncating$type <- rep('pkd2_non_truncating', nrow(pkd2_non_truncating))

#Other genes: self-explanatory 
other_genes <- f[which(f$tier == 'TIER1' & f$genomic_feature_hgnc != 'PKD1'),]
other_genes <- other_genes[which(other_genes$tier == 'TIER1' & other_genes$genomic_feature_hgnc != 'PKD2'),]
other_genes <- other_genes[!duplicated(other_genes$participant_id),]  
other_genes$type <- rep('other_genes', nrow(other_genes))

#NMD: no TIER1 variants or no PKD1/2 variants
nmd <- f[which(f$tier != 'TIER1'),]
nmd <- nmd[which(nmd$genomic_feature_hgnc != 'PKD1'),]
nmd <- nmd[which(nmd$genomic_feature_hgnc != 'PKD2'),]
nmd <- nmd[!duplicated(nmd$participant_id),]  
nmd$type <- rep('nmd', nrow(nmd))

#merge all table type together
all_mutations <- rbind(nmd, other_genes, pkd1_non_truncating, pkd2_non_truncating, pkd1_truncating, pkd2_truncating)

#take highest tier per person for tier type and then by consequence
all_mutations <- subset(all_mutations, as.logical(ave(tier, participant_id, FUN = function(x) x== min(x))))
all_mutations$type <- ordered(all_mutations$type, levels=rev(c("pkd1_truncating", "pkd2_truncating", "pkd1_non_truncating",  "pkd2_non_truncating", "other_genes", "nmd")))
all_mutations <- do.call(rbind, lapply(split(all_mutations, all_mutations$participant_id), function(x) x[which.max(x$type),]))

#merge on participant_id keeping all the values in the age of ESRF table (add all.y if you want to keep everyone)
mutation_age_esrf <- merge(final_hes_combined, all_mutations, by = "participant_id", all.y = T)

#subset table into what is needed
pkd_kaplan <- subset(mutation_age_esrf, select=c(CKD_grade, age_at_last_followup_or_ESRF, type))
#turn CKD grades into 0 or 1
pkd_kaplan$CKD_grade <- ifelse(pkd_kaplan$CKD_grade %in% c(0:4), 1, 5)
pkd_kaplan$CKD_grade[pkd_kaplan$CKD_grade == 5] <- 2
names(pkd_kaplan)[names(pkd_kaplan)=="CKD_grade"] <- "status"
names(pkd_kaplan)[names(pkd_kaplan)=="age_at_last_followup_or_ESRF"] <- "time"
names(pkd_kaplan)[names(pkd_kaplan)=="type"] <- "mutation"
pkd_kaplan <- as.data.frame(pkd_kaplan)
pkd_kaplan$time <- as.numeric(pkd_kaplan$time)
write.table(pkd_kaplan, "~/re_gecip/renal/oalavijeh/projects/phenotypes/stratified/pkd_kaplan_bioinformaticallyv16_svs_as_truncating.txt", col.names = T, row.names = F, quote = F, sep = '\t')
library(survival)
library(ggfortify)
library(survminer)

fit <- survfit(Surv(time, status) ~ mutation, data = pkd_kaplan)
autoplot(fit, conf.int = FALSE, censor = TRUE, main = "renal survival in 100K cystic renal disease cohort (n=419)", xlab ="Age at ESRF", ylab = "survival", surv.geom = 'step')

#see total graph using surviner 
ggsurvplot(fit, data = pkd_kaplan)

#customise graph in survminer 
pp = ggsurvplot(
  fit, 
  data = pkd_kaplan, 
  size = 1,                  # change line size
  risk.table = TRUE,         # Add risk table
  risk.table.col = "strata", # Risk table color by groups
  risk.table.height = 0.31,  # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  ylab = "Renal Survival",
  xlab= "Age/years",
  title = "100K cystic disease cohort", 
)
