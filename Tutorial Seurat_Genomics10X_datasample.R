## Tutorial Seurat: analysis of 10XGenomics sample data 

#load the necessary libraries
library(Seurat)
library(patchwork)
library(dplyr)

#read data
pbmc.data <- Read10X(data.dir = "/Users/carloleonardi/Downloads/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#The [[ operator can add columns to object metadata. This is a great place to stash QC stats. or other measurements stored in the metadata.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot, visualize QC metrics in order then to filter out the outliers
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter cells that have unique feature counts over 2,500 or less than 200 & cells that have >5% mitochondrial counts

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data 
#LogNormalize --> normalize the expression of each cell on total expression, multiply by 10k and logtransform it

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# then we scale the data
# scaling the data is a standard pre-processing step prior to dimensional reduction techniques 

#scale data function --> 1) shifts the expression of each gene, so that mean expr across cells is 0
# 2) scale the expression of each gene, so that variance is 1 --> equal weight, high expr genes do not dominate

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run a PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#look at the PCA in different ways:
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#print list of genes and pc value, easy to discriminate which gene contributes more or less to pc
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# classical pca scatter plot
DimPlot(pbmc, reduction = "pca")

#create a heatmap for the exploration of the primary sources of heterogeneity in a dataset, useful to select which PCs to further analyze
#heatmap is cool, similar info to VizDimLoading but more useful 

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the dimensionality of the data
# to overcome noise --> seurat cluster based on PCA scores (each PC is a metafeature of a correlated feature set)
# how many pc? --> permute a subset of the data and rerun PCA --> construct null distr for features --> strong enrich for low p-value features --> significant pc
# can take some time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

# otherwise: ranking of principle components based on the percentage of variance --> after the drop you have a PC limit
ElbowPlot(pbmc) #you see a dropa after 10, then consider PC 1-10

# clustering the cells
#in brief: we use KNN or SNN (part 1) and then Louivan Algorithm (part 2)

#part 1: cells as a graph structure --> cells are nodes, edges between similar cells in the pc space --> then partition the graph in cliques
# --> knn using euclidean space in the pca space, then weight the local graph based on jaccard similairty btw two nodes (in common genes over all genes)
pbmc <- FindNeighbors(pbmc, dims = 1:10)

#part 2: use the Louvaine algorithm to group cells together --> contains a resolution parameter that sets the ‘granularity’ of the downstream clustering
pbmc <- FindClusters(pbmc, resolution = 0.5)

# use Idents function to see the clusters, here first 5
head(Idents(pbmc), 5)

#Run UMAP, tSNE or other dimensionality-reduction methods

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)

#I found high perplexity helps with messy data
pbmc <- RunTSNE(pbmc, dims = 1:10, perplexity = 80)

#need to install FIt-SNE, look for that

DimPlot(pbmc, reduction = "tsne")

# then, find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#chiedi a giulia, for loop in order to compute all the clusters 
i <- c(1:8)
for i in c(1:length(i))
cluster[i].markers <- FindMarkers(pbmc, ident.1 = i, min.pct = 0.25)

#in alternativa

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)

cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)

cluster6.markers <- FindMarkers(pbmc, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)

#otherwise another way: you can distinguish other markers specifying ident.2, not sure if this can be done in order to substitute the loop
#chiedi a giulia
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#then: test the significance of the markers, using test.use = "roc"
#again: how this can be done in a loop way? 
#ask giulia
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#look the expression of the markers using visualization tools 
#VlnPlot --> expression probability distributions across clusters, FeaturePlot --> visualizes feature expression on a tSNE or PCA plot

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))


FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

#plot raw data as well

VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#another useful tool: create heatmap for given cell and features
#here selecting for top 10 markers

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#in this case: clusters are known (markers are gene signatures, we can recognise a cell pop)
#let's name them and plot in the umap

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#don't forget to save 
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

#further questions for giulia:
# reading this: https://www.embopress.org/doi/10.15252/msb.20188746 
#what about compositional analysis? and trajectory analysis (and in definying trajecotry, what about tracing metastable states in pseudotime?)

i <- c(1:8)
for (i in c(1:length(i))))
eval(parse(paste("cluster[",i,"].markers" <- FindMarkers(pbmc, ident.1 = ",i,", min.pct = 0.25),sep="")))



