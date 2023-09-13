#############################################################################################################################

#USE Bedtools to filter individual SV vcf files #################################################################

# This script subsets VCFs by coordinate and prints out sample ID and genotype

# INPUT:
        # 1: Tab-delimited list of chromosomal coordinates (regions.txt)
        # 2. A list of full-paths to VCF files of interest (vcf_list.txt)

# OUTPUT:
        # 1: A tab-delimted file of sample ID, variant information and quality (results.txt)

# STEPS:
        # 1: tabix to subset each VCF by the coordinates given
        # 2: bcftools to split multi-allelic variants across multiple lines
        # 3: bcftools to filter for PASS variants and variants with a coverage >10 and genotype quality >15
        # 4: bcftools to print out the sample ID, and variant and quality information as a text file

# IMPORTANT!
        # 1: Make sure the chromosome nomenculature is correct for GRCh37 and GRCh38
        # 2: Do not mix genome builds in VCF file list
        # 3: Change the file paths to suit your working environment
	# 4: This script will only extract SNVs, if you need indels please comment out the '-i 'MIN(FMT/DP)>10 & MIN(FMT/GQ)>15' flag
	# 5: This script will only be function with Illumina germine variant calls - please use the relevent script to process somatic variants

#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o sv_%J.stdout
#BSUB -e sv_%J.stderr
#BSUB -cwd /re_gecip/renal/oalavijeh/projects/pkd_sv/data/exon_wide/pkd
#BSUB -R "rusage[mem=10000] span[hosts=1]"

module load bio/BEDTools/2.27.1-foss-2018b
module load bio/BCFtools/1.11-GCC-8.3.0

COHORT=pkd
GENE=exonwide

VCFS="/re_gecip/renal/oalavijeh/projects/pkd_sv/data/exon_wide/pkd/pkd_vcf_list.txt"
REGIONS="/renal/oalavijeh/projects/pkd_sv/data/workflow/pkd_genecode.bed"
EXONS="/re_gecip/renal/mchan/gencode_v29_exons.bed"
OUTPATH="/re_gecip/renal/oalavijeh/projects/pkd_sv/data/exon_wide/pkd"


while read -r vcf; do
	vfile="${vcf##*/}"
        OUTFILE="${OUTPATH}/${vfile}.${COHORT}_${GENE}.vcf"
	bedtools intersect -a ${vcf} -b ${EXONS} -u -header | \
	bcftools view -f PASS -o ${OUTFILE}
done < ${VCFS}

cd ${OUTPATH}
ls *.${COHORT}_${GENE}.vcf > ${COHORT}_${GENE}_vcf.txt
