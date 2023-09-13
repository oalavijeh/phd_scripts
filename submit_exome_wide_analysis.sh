#! /usr/bin/env bash

#BSUB -q long
#BSUB -P re_gecip_renal
#BSUB -o sv_exonwide_%J.stdout
#BSUB -e sv_exonwide_%J.stderr
#BSUB -cwd /re_gecip/renal/oalavijeh/projects/pkd_sv/analysis/exomewide_nmd
#BSUB -R "rusage[mem=10000] span[hosts=1]"
#BSUB -W 360:0

module load toolchain/foss/2019b
module load lang/R/3.6.2-foss-2019b
module load bio/HTSlib/1.9-GCC-8.3.0

Rscript /re_gecip/renal/oalavijeh/projects/pkd_sv/scripts/new_workflowv2_2022/SV_exome_wide_analysis_maf0.001_nmd.R
