library( infercnv )
library( Matrix )

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript infercnv_analysis.R list_patients_folders.txt\n")
  quit()
}

wd_list <- readLines(args[1])

for(wd in unique(wd_list)){
  
  setwd(wd)
  
  cellranger_outs_folder_negative <- NA
  cell.annotation.name <- paste0("cellAnnotations_filtered_",basename(wd),".txt")
  counts.sparse.matrix.name <- paste0("sc.10x.counts_",basename(wd),".RData")
  
  out.folder.name <- file.path(wd,"results")
  
  load(counts.sparse.matrix.name)
  
  # create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat.matrix.sparse,
                                      annotations_file=cell.annotation.name,
                                      delim="\t",
                                      gene_order_file="gene_ordering_file.txt",
                                      ref_group_names=c("normal"))
  #ref_group_names=NULL)
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=out.folder.name,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,
                               HMM=F
  )

}
