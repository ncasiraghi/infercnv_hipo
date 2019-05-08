#!/bin/sh 
#BSUB -J infercnv
#BSUB -q night 
#BSUB -e infercnv.bsub.log 
#BSUB -o infercnv.bsub.txt 

module load R/3.5.1

Rscript infercnv_analysis_bsub.R
