library( Seurat )

patients <- list.files("/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/MM_scRNAseq_data/K43R/Patients",full.names = T)

clean_barcodes <- list()

for(p in unique(patients)){
  message(basename(p))
  
  tumor_files <- list.files(p,pattern = paste0("Tumor_"),full.names = T)
  
  get_barcodes <- function(tumor_files, pattern){
    rds <- grep(tumor_files,pattern = pattern,value = T)
    if(length(rds)>0){
      tumor <- readRDS(rds)
      tumorDF <- tumor@meta.data[grep(rownames(tumor@meta.data),pattern = "_T[[:digit:]]_"),]
      return(rownames(tumorDF))
    }
  }
  
  pre <- get_barcodes(tumor_files = tumor_files,pattern = "_pre.rds")
  post <- get_barcodes(tumor_files = tumor_files,pattern = "_post.rds")

  get_clean_barcode <- function(bc){
    sample <- unlist(strsplit(bc,split = '_'))[2]
    barcode <- unlist(strsplit(bc,split = '_'))[3]
    return(paste0(barcode,'-1pos_',sample))
  }
  
  clean_barcodes[[basename(p)]] <- unique(unlist(lapply(unique(c(pre,post)),FUN = get_clean_barcode)))
}



