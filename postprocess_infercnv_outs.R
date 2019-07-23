library( phylogram )
library( dendextend ) 

infercnv.observations <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/infercnv_outs/clean_barcodes/TS6ZX9/results/infercnv.preliminary.observations_dendrogram.txt"

# get and visualize the dendrogram

dend <- as.hclust(read.dendrogram(infercnv.observations))
plot(dend,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

# cut the dendrogram 
NCLONES <- 4

tp <- color_branches(dend, k = NCLONES)
plot(tp,leaflab = 'none')

g <- cutree(dend,k = NCLONES)

# adjust barcodes
clean_barcode <- function(b){
  id <- unlist(strsplit(b,split = '_'))[2]
  new <- paste0('TS6ZX9_',id,'_',gsub(b,pattern = "-1pos_T1|-1pos_N1",replacement = ""))
  return(new)
}

barcodes <- unlist(lapply(names(g),clean_barcode))

groups <- data.frame(BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

save(groups,file = 'TS6ZX9_Barcodes_Clones.RData',compress = TRUE)
