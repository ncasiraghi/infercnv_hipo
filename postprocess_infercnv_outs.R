library( phylogram )
library( dendextend ) 
library( gplots )
library( data.table )
library( infercnv )

# functions

clean_barcode <- function(b){
  id <- unlist(strsplit(b,split = '_'))[2]
  new <- paste0('TS6ZX9_',id,'_',gsub(b,pattern = "-1pos_T1|-1pos_N1",replacement = ""))
  return(new)
}

## standard output from InferCNV

infercnv.observations <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/N5CC3E/results/infercnv.observations_dendrogram.txt"

## different input

setwd('/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/N5CC3E/results')

infercnv_obj = readRDS('run.final.infercnv_obj')

tumor_expr_data <- infercnv_obj@observation_grouped_cell_indices

bc <- read.delim('../cellAnnotations_filtered_N5CC3E.txt',header = F)
barcodes <- unique(grep(bc$V1,pattern = '-1neg',invert = TRUE,value = T))

tumor_expr_data <- tumor_expr_data[,barcodes]
dend <- hclust(d = dist(t(tumor_expr_data)), method = "ward.D2")

# get and visualize the dendrogram
dend <- as.hclust(read.dendrogram(infercnv.observations))
plot(dend,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

length(dend$labels)

# cut the dendrogram 
NCLONES <- 4

tp <- color_branches(dend, k = NCLONES,groupLabels = T)
plot(tp,leaflab = 'none',horiz = F) # set horiz = TRUE to plot dendrogram as in infercnv heatmap

g <- cutree(dend,k = NCLONES)
barcodes <- unlist(lapply(names(g),clean_barcode))

groups <- data.frame(BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

#save(groups,file = 'TS6ZX9_Barcodes_Clones.RData',compress = TRUE)


## custom analysis from InferCNV final matrix

# Get the the tumor cell matrix data values
input <- '/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/TS6ZX9/results/infercnv.observations.txt'
finalmat <- read.table(file = input,header = T,row.names = 1)
dim(finalmat)

genes <- read.delim('/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/TS6ZX9/gene_ordering_file.txt',as.is = T,stringsAsFactors = F,header = F)
colnames(genes) <- c("gene_id","chr","start","end")

keep_genes <- which(rownames(finalmat) %in% unique(genes$gene_id[which(genes$chr %in% c('chr5','chr6','chr10','chr18'))]))

m <- t(as.matrix(finalmat[keep_genes,]))

# which(apply(m, 1, var)==0)
# m <- m[ , apply(m, 1, var) != 0]
# pca <- prcomp(t(m),scale. = T)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 10)

heatmap.2(m,
          col = my_palette,
          Rowv = TRUE,
          Colv = FALSE,
          trace="none",
          scale="none",
          dendrogram = "row",
          labCol=FALSE,labRow = FALSE)

hc.rows<- hclust(dist(m))

plot(hc.rows,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

# cut the dendrogram 
NCLONES <- 5

tp <- color_branches(hc.rows, k = NCLONES,groupLabels = T)
plot(tp,leaflab = 'none',horiz = F) # set horiz = TRUE to plot dendrogram as in infercnv heatmap

g <- cutree(hc.rows,k = NCLONES)
barcodes <- unlist(lapply(names(g),clean_barcode))

groups <- data.frame(BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

#save(groups,file = 'TS6ZX9_Barcodes_Clones_custom_clustering.RData',compress = TRUE)
