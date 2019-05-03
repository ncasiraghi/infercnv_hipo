library( infercnv )
library( Matrix )

wd <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/"
setwd(wd)

counts.matrix.name <- file.path(wd,"sc.10x.counts_K43R-8YGUU8-T3Z3.matrix")
cell.annotation.name <- file.path(wd,"cellAnnotations_K43R-8YGUU8-T3Z3.txt")
out.folder.name <- file.path(wd,"K43R-8YGUU8-T3Z3")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts.matrix.name,
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
