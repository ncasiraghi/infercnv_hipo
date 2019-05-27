library( Matrix )

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  quit()
}
source(args[1])

setwd(wd)

# Accomodating 10X Data

getCountMat <- function(cellranger_outs_folder){
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
  return( mat )
}

getBarcodeNames <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  return( as.character(barcode.names[,1]) )
}

# get count matrix

mat.pos.list <- lapply( cellranger_outs_folder_positive, getCountMat )
mat.neg.list <- lapply( cellranger_outs_folder_negative, getCountMat )

genes <- unique(unlist(lapply(mat.pos.list, rownames)))

filter_genes <- function(x, genes){
  x <- x[genes,]
  return(x)
}

mat.pos.list <- lapply(mat.pos.list, filter_genes, genes)

# add label for cells from the same sample
for(i in seq(length(mat.pos.list))){
  colnames(mat.pos.list[[i]]) <- paste0(colnames(mat.pos.list[[i]]), paste0("pos_",labels.pos[i]))
}

for(i in seq(length(mat.neg.list))){
  colnames(mat.neg.list[[i]]) <- paste0(colnames(mat.neg.list[[i]]), "neg")
}

check_and_merge <- function(sparse_matrix_list){
  check.cols <- c()
  if(length(sparse_matrix_list) < 2){
    return(sparse_matrix_list[[1]])
  }
  combtab <- combn(x = seq(length(sparse_matrix_list)),m = 2)
  for(i in seq(ncol(combtab))){
    check.cols <- identical(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]]))
  }
  if(all(check.cols)){
    mat.pos <- do.call(cbind, sparse_matrix_list)
    return(mat.pos)
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
save(mat.matrix.sparse, file = counts.sparse.matrix.name,compress = T)

# cell annotations

cellAnnotations <- data.frame(cell_ids=colnames(mat.matrix.sparse),
                              annotation=c(rep(labels.pos,as.numeric(lapply(mat.pos.list,ncol))),
                                           rep("normal",as.numeric(lapply(mat.neg.list,ncol)))
                                           ),
                              stringsAsFactors = F)

# only when mixed positive and negative fractions
# cellAnnotations <- data.frame(cell_ids=c(cells.pos, cells.neg),
#                               annotation=c( rep("T1",length(cells.pos)), rep("N1",length(cells.neg))),
#                               stringsAsFactors = F)

write.table(cellAnnotations, file=cell.annotation.name, col.names = F, row.names = F,quote=F, sep="\t")                   

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
  write.table(gene_ordering_file[,c(5,1:3)], file='gene_ordering_file.txt', col.names = F, row.names = F,quote=F, sep="\t")
}
