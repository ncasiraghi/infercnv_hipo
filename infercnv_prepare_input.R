library( Matrix )

wd <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/"
setwd(wd)

## K43R-ZPMZFJ
# # inputs 
# cellranger_outs_folder_positive <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-ZPMZFJ-T1/outs/"
# cellranger_outs_folder_negative <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-ZPMZFJ-N1/outs/"
# # outputs
# cell.annotation.name <- "cellAnnotations_K43R-ZPMZFJ-T1N1.txt"
# #cell.annotation.name <- "cellAnnotations_K43R-ZPMZFJ-T1N1.mixed.txt"
# counts.matrix.name <- "sc.10x.counts_K43R-ZPMZFJ-T1N1.matrix"

## K43R-8YGUU8 T2Z2
# # inputs
# cellranger_outs_folder_positive <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-T2/outs/"
# cellranger_outs_folder_negative <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-Z2/outs/"
# # outputs
# cell.annotation.name <- "cellAnnotations_K43R-8YGUU8-T2Z2.txt"
# counts.matrix.name <- "sc.10x.counts_K43R-8YGUU8-T2Z2.matrix"

## K43R-8YGUU8 T3Z3
# # inputs
# cellranger_outs_folder_positive <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-T3/outs/"
# cellranger_outs_folder_negative <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-Z3/outs/"
# # outputs
# cell.annotation.name <- "cellAnnotations_K43R-8YGUU8-T3Z3.txt"
# counts.matrix.name <- "sc.10x.counts_K43R-8YGUU8-T3Z3.matrix"

## K43R-8YGUU8 T4Z4
# inputs
cellranger_outs_folder_positive <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-T4/outs/"
cellranger_outs_folder_negative <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3/K43R-8YGUU8-Z4/outs/"
# outputs
cell.annotation.name <- "cellAnnotations_K43R-8YGUU8-T4Z4.txt"
counts.sparse.matrix.name <- "sc.10x.counts_K43R-8YGUU8-T4Z4.RData"


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
mat.pos <- getCountMat(cellranger_outs_folder = cellranger_outs_folder_positive)
mat.neg <- getCountMat(cellranger_outs_folder = cellranger_outs_folder_negative)

# add label for cells
colnames(mat.pos) <- paste0(colnames(mat.pos),"pos")
colnames(mat.neg) <- paste0(colnames(mat.neg),"neg")

if(identical(rownames(mat.pos),rownames(mat.neg))){
  mat.matrix.sparse <- cbind(mat.pos, mat.neg)
} else {
  message("mat.pos and mat.neg have different number or rownames")
  stop()
}

# < check >
dim(mat.pos)
dim(mat.neg)
dim(mat.matrix.sparse)

# merge pos and neg
save(mat.matrix.sparse, file = counts.sparse.matrix.name,compress = T)

# cell annotations
cells.pos <- paste0(getBarcodeNames(cellranger_outs_folder = cellranger_outs_folder_positive),"pos")
cells.neg <- paste0(getBarcodeNames(cellranger_outs_folder = cellranger_outs_folder_negative),"neg")

cellAnnotations <- data.frame(cell_ids=c(cells.pos, cells.neg),
                              annotation=c( rep("tumor",length(cells.pos)), rep("normal",length(cells.neg))),
                              stringsAsFactors = F)

# cellAnnotations <- data.frame(cell_ids=c(cells.pos, cells.neg),
#                               annotation=c( rep("T1",length(cells.pos)), rep("N1",length(cells.neg))),
#                               stringsAsFactors = F)



write.table(cellAnnotations[,1:2], file=cell.annotation.name, col.names = F, row.names = F,quote=F, sep="\t")                   

# genes ordering
gene_ordering_file <- read.delim(file = "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed",stringsAsFactors = F,skip = 1,header = F)
gene_ordering_file <- gene_ordering_file[which(gene_ordering_file$V1 != "MT"),]
gene_ordering_file$index <- gene_ordering_file$V1

gene_ordering_file$index[which(gene_ordering_file$index == "X")] <- 23
gene_ordering_file$index[which(gene_ordering_file$index == "Y")] <- 24
gene_ordering_file$index <- as.numeric(gene_ordering_file$index)
gene_ordering_file$V1 <- paste0("chr",gene_ordering_file$V1)
gene_ordering_file <- gene_ordering_file[with(gene_ordering_file,order(index)),]
write.table(gene_ordering_file[,c(5,1:3)], file='gene_ordering_file.txt', col.names = F, row.names = F,quote=F, sep="\t")
