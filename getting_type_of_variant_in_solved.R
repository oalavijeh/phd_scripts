# script to ascertain type of variant listed as solved 
# load labkey
library(Rlabkey)

# Set the baseURL

labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")

# Write your SQL query here - if a new metric is in a non-listed table in the FROM section - needs defining
# SELECT are coloumns of interest.
# WHERE are the search function - disease, type of file (Array GenotypingGenomic/Repeat/Standard/Structural VCF or BAM),build(GRCh38 or GRCh37),participant (Proband, Relative),pick europeans with >0.9

query <- "SELECT rd.participant_id, rd.participant_type, par.consanguinity, par.participant_phenotypic_sex, par.year_of_birth, par.mother_affected, par.father_affected, par.full_brothers_affected, par.full_sisters_affected, rd.normalised_specific_disease, rd.plate_key, gc.acmg_classification, gc.case_solved_family, gc.additional_comments, gc.gene_name, gc.assembly, gc.chromosome, gc.position, gc.reference, gc.alternate
FROM rare_disease_analysis AS rd
LEFT JOIN gmc_exit_questionnaire gc
ON rd.participant_id = gc.participant_id
LEFT JOIN participant AS par
ON gc.participant_id = par.participant_id
WHERE rd.normalised_specific_disease = 'Cystic kidney disease'
OR rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
OR rd.normalised_specific_disease = 'Proteinuric renal disease'
OR rd.normalised_specific_disease = 'Unexplained kidney failure in young people'
OR rd.normalised_specific_disease = 'Extreme early-onset hypertension'
OR rd.normalised_specific_disease = 'Familial IgA nephropathy and IgA vasculitis'
OR rd.normalised_specific_disease = 'Renal tract calcification (or Nephrolithiasis or nephrocalcinosis)'
OR rd.normalised_specific_disease = 'Renal tubular acidosis'
OR rd.normalised_specific_disease = 'Atypical haemolytic uraemic syndrome'
OR rd.normalised_specific_disease = 'Primary membranoproliferative glomerulonephritis'
OR rd.normalised_specific_disease = 'Familial haematuria'
AND rd.participant_type = 'Proband'"
         

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v16_2022-10-13",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice



)

y <-data.frame(mysql)
y<- y[!duplicated(y$participant_id),]
pkd <-subset(y, participant_type =='Proband' & normalised_specific_disease == "Cystic kidney disease")
#pkd[is.na(pkd)] <- "Unknown"
solved_pkd <- subset(pkd, case_solved_family == "yes")
unsolved_pkd <- subset(pkd, case_solved_family == "no")
#get tiering list
# Set the baseURL
#Write your SQL query here - if a new metric is in a non-listed table in the FROM section - needs defining
# SELECT are coloumns of interest.
# WHERE are the search function - disease, type of file (Array GenotypingGenomic/Repeat/Standard/Structural VCF or BAM),build(GRCh38 or GRCh37),participant (Proband, Relative),pick europeans with >0.9

query <- "SELECT td.participant_id, td.sample_id, td.tier,td.assembly, td.chromosome,td.position,td.reference,td.alternate,td.genomic_feature_hgnc, td.consequence_type
          FROM tiering_data as td
          WHERE td.phenotype = 'Cystic kidney disease'
          AND td.participant_type = 'Proband'"

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v16_2022-10-13",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice
  
  
  
)

#match tier list to solved list to get consequence
pkd1or2_tier<-subset(mysql, mysql$genomic_feature_hgnc == "PKD1" | 
                       mysql$genomic_feature_hgnc == "PKD2" |
                       mysql$genomic_feature_hgnc == "DNAJB11"|
                       mysql$genomic_feature_hgnc == "GANAB"|
                       mysql$genomic_feature_hgnc == "BBS1"|
                       mysql$genomic_feature_hgnc == "HNF1B"|
                       mysql$genomic_feature_hgnc == "COL4A4"|
                       mysql$genomic_feature_hgnc == "OFD1"|
                       mysql$genomic_feature_hgnc == "TMEM67"|
                       mysql$genomic_feature_hgnc == "UMOD"|
                       mysql$genomic_feature_hgnc == "WT1"|
                       mysql$genomic_feature_hgnc == "SDCCAG8"|
                       mysql$genomic_feature_hgnc == "SALL1"|
                       mysql$genomic_feature_hgnc == "FAN1")
                       
solved_by_type <- merge(solved_pkd, pkd1or2_tier, all.x = T, all.y = F)
#split into PKD1 and then PKD2 
#create hierarchical preference 
i <- order(match(solved_by_type$consequence_type, c("stop_gained", "stop_lost", 
                                                    "splice_donor_variant", 
                                                    "splice_acceptor_variant", 
                                                    "frameshift_variant","inframe_deletion", 
                                                    "inframe_insertion", "missense_variant", 
                                                    "splice_region_variant", "synonymous_variant", 
                                                    "5_prime_UTR_variant", "3_prime_UTR_variant",
                                                    "non_coding_transcript_exon_variant", 
                                                    "intron_variant", "NMD_transcript_variant",
                                                    "non_coding_transcript_variant", 
                                                    "upstream_gene_variant", "downstream_gene_variant",
                                                    "2KB_downstream_variant", "2KB_upstream_variant",
                                                    "NA")))
v2solved_by_type <- solved_by_type[i[!duplicated(solved_by_type$participant_id[i])],]
#sort gene_name and get table
replace_list <- c("TSC2;PKD1.*" = "TSC2/PKD1","PKD1;.*"= "PKD1","PKD2;.*" = "PKD2",
                  "DNAJB11.*" = "DNAJB11","PKHD1.*" = "PKHD1","BBS1.*"="BBS1")
for (item in names(replace_list)) {
 v2solved_by_type$gene_name <- gsub(item, replace_list[item], v2solved_by_type$gene_name)
}
#get table of solved 
solved_cons_gene <- as.data.frame(table(v2solved_by_type$gene_name, v2solved_by_type$consequence_type))
solved_cons_gene <- subset(solved_cons_gene, solved_cons_gene$Freq >0)

write.table(v2solved_by_type, "~/re_gecip/renal/oalavijeh/projects/phenotypes/final_phenotypes/solved_with_variant_consequence_v16.txt", col.names = T, row.names = F, sep = '\t', quote = F)
