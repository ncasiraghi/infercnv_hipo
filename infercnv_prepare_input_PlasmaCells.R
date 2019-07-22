library( Matrix )
library( Seurat )

# load the plasma cell file
plasma_cells_file <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/MM_scRNAseq_data/K43R/ME_only/Plasma_cells.rds"
pc <- readRDS(plasma_cells_file)

# select patients from wich extract plasma cells
patients_to_keep <- unique(pc@meta.data$patient)
patients <- list.files("/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/MM_scRNAseq_data/K43R/Patients",full.names = T)
patients <- patients[which(basename(patients) %in% patients_to_keep)]

# outs
main_outdir <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/PlasmaCells/"
cell.annotation.name <- paste0("cellAnnotations_filtered_PlasmaCells.txt")
counts.sparse.matrix.name <- paste0("sc.10x.counts_PlasmaCells.RData")

#wd <- patients[1]
keepBarcodes <- list()

for(wd in patients){
  
  message(wd)
  
  cellranger_all_results <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3"
  
  cellranger_outs_folder_positive <- file.path(list.files(cellranger_all_results, pattern = basename(wd),full.names = T),'outs/')
  
  getlabel <- function(x,pattern=basename(wd)){
    return(unlist(strsplit(grep(unlist(strsplit(x,split = '/')),pattern = pattern,value = T),split = paste0(pattern,'-')))[2])
  }
  
  labels.pos <- unlist(lapply(cellranger_outs_folder_positive, getlabel))
  
  allBarcodes <- grep(unique(rownames(pc@meta.data)),pattern = basename(wd),value = T)
  
  for(lab in labels.pos){
    bc <- allBarcodes[grep(allBarcodes,pattern = lab)]
    
    clean <- function(x,TAG){
      return( paste0(unlist(strsplit(x,split = paste0(TAG,"_")))[2],'-1') )
    }
    
    keepBarcodes[[paste0(basename(wd),"-",lab)]] <- unlist(lapply(bc, clean, TAG=lab))
    
  }
}

# Accomodating 10X Data

getCountMat <- function(cellranger_outs_folder, barcodes_to_keep){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix_dir, "features.tsv.gz")
  matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  key <- gsub(basename(gsub(cellranger_outs_folder,pattern = '/outs/',replacement = '')),pattern = 'K43R-',replacement = '')
  
  if(!is.null(barcodes_to_keep[[key]])){
    mat <- mat[,barcodes_to_keep[[key]],drop=FALSE]
    
    return( mat )
    
  } else {
    return( NULL )
  }
}

# get count matrix
cellranger_all_results <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3"

cellranger_outs_folder_negative <- NA
cellranger_outs_folder_positive <- file.path(list.files(cellranger_all_results,full.names = T),'outs/')

mat.pos.list_raw <- lapply( cellranger_outs_folder_positive, getCountMat, barcodes_to_keep=keepBarcodes)
mat.pos.list <- plyr::compact(mat.pos.list_raw)

if(is.na(cellranger_outs_folder_negative)){
  load("/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/HumanCellAtlas_plasmacells/HCA_PCs_counts_ensmbl.RData")
  mat.neg.list <- list()
  mat.neg.list[[1]] <- hca
} else {
  mat.neg.list <- lapply( cellranger_outs_folder_negative, getCountMat )
}

# add label for cells from the same sample
for(i in seq(length(mat.pos.list))){
  colnames(mat.pos.list[[i]]) <- paste0(colnames(mat.pos.list[[i]]),"pos")
}

for(i in seq(length(mat.neg.list))){
  colnames(mat.neg.list[[i]]) <- paste0(colnames(mat.neg.list[[i]]), "neg")
}

check_and_merge <- function(sparse_matrix_list){
  if(length(sparse_matrix_list) < 2){
    return(sparse_matrix_list[[1]])
  }
  
  l <- list()
  for(i in seq_len(length(sparse_matrix_list))){
    l[[i]] <- rownames(sparse_matrix_list[[i]])
  }
  red <- Reduce(intersect,l)
  
  filter_genes <- function(x, genes){
    x <- x[genes,,drop=FALSE]
    return(x)
  }
  
  sparse_matrix_list <- lapply(sparse_matrix_list, filter_genes, genes=red)
  
  check <- list()
  for(i in seq_len(length(sparse_matrix_list))){
    check[[i]] <- rownames(sparse_matrix_list[[i]])
  }
  
  check <- do.call(cbind,check)
  
  if(length(apply(check,MARGIN = 1,FUN = unique)) == length(red)){
    outmat <- do.call(cbind, sparse_matrix_list)
    return(outmat)
  } else {
    return(NULL)
  }
  
}

mat.pos <- check_and_merge(sparse_matrix_list = mat.pos.list)
mat.neg <- check_and_merge(sparse_matrix_list = mat.neg.list)

mat.matrix.sparse <- check_and_merge(sparse_matrix_list = list(mat.pos,mat.neg))

dim(mat.pos)
dim(mat.neg)
dim(mat.matrix.sparse)

# merge pos and neg
for(i in seq_len(ncol(mat.matrix.sparse))){
  colnames(mat.matrix.sparse)[i] <- paste0(colnames(mat.matrix.sparse)[i],"-",as.character(i))
}

save(mat.matrix.sparse, file = file.path(main_outdir,counts.sparse.matrix.name),compress = T)

# cell annotations

cellAnnotations <- data.frame(cell_ids=colnames(mat.matrix.sparse),
                              annotation=c(rep("PlasmaCells",ncol(mat.pos)),
                                           rep("normal",ncol(mat.neg))
                              ),
                              stringsAsFactors = F)

# only when mixed positive and negative fractions
# cellAnnotations <- data.frame(cell_ids=c(cells.pos, cells.neg),
#                               annotation=c( rep("T1",length(cells.pos)), rep("N1",length(cells.neg))),
#                               stringsAsFactors = F)

write.table(cellAnnotations, file=file.path(main_outdir,cell.annotation.name), col.names = F, row.names = F,quote=F, sep="\t")                   

# genes ordering
if( TRUE ){
  gene_ordering_file <- read.delim(file = "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed",stringsAsFactors = F,skip = 1,header = F)
  gene_ordering_file <- gene_ordering_file[which(gene_ordering_file$V1 != "MT"),]
  gene_ordering_file$index <- gene_ordering_file$V1
  
  gene_ordering_file$index[which(gene_ordering_file$index == "X")] <- 23
  gene_ordering_file$index[which(gene_ordering_file$index == "Y")] <- 24
  gene_ordering_file$index <- as.numeric(gene_ordering_file$index)
  gene_ordering_file$V1 <- paste0("chr",gene_ordering_file$V1)
  gene_ordering_file <- gene_ordering_file[with(gene_ordering_file,order(index)),]
  write.table(gene_ordering_file[,c(5,1:3)], file=file.path(main_outdir,'gene_ordering_file.txt'), col.names = F, row.names = F,quote=F, sep="\t")
}
