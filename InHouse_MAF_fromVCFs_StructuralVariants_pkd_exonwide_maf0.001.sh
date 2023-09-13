#!/usr/bin/env bash

#BSUB -q long
#BSUB -P re_gecip_renal
#BSUB -o svmaf_%J.stdout
#BSUB -e svmaf_%J.stderr
#BSUB -cwd /re_gecip/renal/oalavijeh/projects/pkd_sv/analysis/exomewide_nmd
#BSUB -R "rusage[mem=10000] span[hosts=1]"

COHORT=maf0.001_nmd
GENE=exonwide

echo $'\n'"["`date`"]: Job started."

##load modules
module load lang/Perl/5.30.0-GCCcore-8.3.0

SCRIPT_PATH="/re_gecip/renal/oalavijeh/projects/pkd_sv/scripts/new_workflow"
GENESFILE="/re_gecip/renal/mchan/gencode_v29_exons.bed"
SAMPLE_PATH="/re_gecip/renal/oalavijeh/projects/pkd_sv/data/exon_wide/pkd"
VCFSFILE="/re_gecip/renal/oalavijeh/projects/pkd_sv/data/exon_wide/pkd/maf0.001/${COHORT}_${GENE}_vcf.txt"
ID="${COHORT}_${GENE}"

##Run perl script
perl "${SCRIPT_PATH}/InHouse_MAF_fromVCFs_StructuralVariants.pl" ${SAMPLE_PATH} ${VCFSFILE} ${ID} ${GENESFILE}

echo $'\n'"["`date`"]: InHouse MAF script is complete!!"
