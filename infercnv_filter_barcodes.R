library( Seurat )

patients <- list.files("/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/MM_scRNAseq_data/K43R/Patients",full.names = T)

for(p in unique(patients)){
  message(basename(p))
  
  
  tumor_file <- list.files(p,pattern = paste0("Tumor_",basename(p),".rds"),full.names = T)
  tumor <- readRDS(tumor_file)
  
  cat(unique(tumor@meta.data),'\n')
}