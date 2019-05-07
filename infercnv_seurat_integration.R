library( dplyr )
library( Seurat )

wd <- "/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/InferCNV/"
setwd(wd)

# sample.name <- "K43R-ZPMZFJ-T1"
# infercnv_folder_outout <- "K43R-ZPMZFJ-T1N1"

sample.name <- "K43R-8YGUU8-T2"
infercnv_folder_outout <- "K43R-8YGUU8-T2Z2"

# sample.name <- "K43R-8YGUU8-T3"
# infercnv_folder_outout <- "K43R-8YGUU8-T3Z3"

# sample.name <- "K43R-8YGUU8-T4"
# infercnv_folder_outout <- "K43R-8YGUU8-T4Z4"

cellranger_outs_folder_positive <- file.path("/icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3",sample.name,"/outs/filtered_feature_bc_matrix/")
seu.data <- Read10X(cellranger_outs_folder_positive)

seu <- CreateSeuratObject(counts = seu.data, project = sample.name, min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")

VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object = seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#seu <- subset(x = seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
seu <- subset(x = seu, subset = percent.mt < 10)

seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)

# scaling data
all.genes <- rownames(x = seu)
seu <- ScaleData(object = seu, features = all.genes)

seu <- RunPCA(object = seu, features = VariableFeatures(object = seu))

print(x = seu[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object = seu, dims = 1:2, reduction = "pca")

DimPlot(object = seu, reduction = "pca")

DimHeatmap(object = seu, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = seu, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(object = seu)

seu <- FindNeighbors(object = seu, dims = 1:10)
seu <- FindClusters(object = seu, resolution = 0.5)
head(x = Idents(object = seu), 5)

seu <- RunUMAP(object = seu, dims = 1:10)
seu <- RunTSNE(object = seu, dims = 1:10)

DimPlot(object = seu, reduction = "umap")
DimPlot(object = seu, reduction = "tsne")

save(seu, file = paste0(sample.name,".RData"),compress = T)

# integrate with infercnv output

# load(paste0(sample.name,".RData"))

head(seu@reductions$umap@cell.embeddings)
head(seu@reductions$tsne@cell.embeddings)
head(seu@meta.data)

FeaturePlot(object = seu,reduction = "tsne",features = "CCND1")

library( phylogram )

dend <- as.hclust( read.dendrogram(file.path(infercnv_folder_outout,"infercnv.observations_dendrogram.txt") ))
class(dend)

plot(dend,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance",main = sample.name)

# cut the tree 
g <- cutree(dend,k = 4)

groups <- data.frame(barcode=gsub(names(g),pattern = "-1pos",replacement = ""),
                     cluster=as.numeric(g),
                     stringsAsFactors = F)

# tsne
tsnetab <- as.data.frame(seu@reductions$tsne@cell.embeddings,stringsAsFactors = F)
tsnetab$barcode <- rownames(tsnetab)
head(tsnetab)
a <- merge(x = tsnetab,y = groups,by = "barcode",all = F)

a$col <- adjustcolor(a$cluster,alpha.f = .6)

par(pty="s",mfrow=c(1,2))
#plot(x = a$tSNE_1,y = a$tSNE_2,pch=16,col=grey(.3,.3))
plot(x = a$tSNE_1,y = a$tSNE_2,pch=20,col=a$col)

library( ape )

colors <- unique(a[,c("cluster","col")])
colors <- colors$col[with(colors, order(cluster))]
plot(as.phylo(dend), tip.color = colors[g],label.offset = 1, cex = 0.1,direction = "rightwards")
axisPhylo(side = 2)

# umap
umaptab <- as.data.frame(seu@reductions$umap@cell.embeddings,stringsAsFactors = F)
umaptab$barcode <- rownames(umaptab)
head(umaptab)
a <- merge(x = umaptab,y = groups,by = "barcode",all.x = T)

par(pty="s",mfrow=c(1,2))
plot(x = a$UMAP_1,y = a$UMAP_2,pch=16,col=grey(.3,.3))
plot(x = a$UMAP_1,y = a$UMAP_2,pch=16,col=adjustcolor(a$cluster,alpha.f = .3))



