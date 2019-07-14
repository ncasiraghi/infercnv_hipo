#!/bin/sh 
#BSUB -J infercnv
#BSUB -q verylong 
#BSUB -e bsub_infercnv.log 
#BSUB -o bsub_infercnv.txt 

module load R/3.5.1

Rscript /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_hipo/infercnv_analysis.R /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/ZPMZFJ/config_infercnv.R
