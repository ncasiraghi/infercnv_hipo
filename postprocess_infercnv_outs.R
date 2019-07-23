library( phylogram )
library( dendextend ) 

infercnv_folder_outout <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/TS6ZX9/results"

# get and visualize the dendrogram

dend <- as.hclust( read.dendrogram(file.path(infercnv_folder_outout,"infercnv.observations_dendrogram.txt") ))
class(dend)

plot(dend,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

# cut the dendrogram 

NCLONES <- 4

tp <- color_branches(dend, k = NCLONES)
plot(tp,leaflab = 'none')

g <- cutree(dend,k = NCLONES)
groups <- data.frame(barcode=gsub(names(g),pattern = "pos_N1|pos_T1",replacement = ""),
                     cluster=as.numeric(g),
                     stringsAsFactors = F)
