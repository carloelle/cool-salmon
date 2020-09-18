
##### new assignment #####
library(Seurat)
library(patchwork)
library(dplyr)

pbmc.data <- Read10X(data.dir = "/Users/carloleonardi/Downloads/SCdata_share_Escobar")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 100)

#The [[ operator can add columns to object metadata. This is a great place to stash QC stats. or other measurements stored in the metadata.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

s.genes<-readLines('/Users/carloleonardi/Desktop/R/cool-salmon/cell_cycle_genes_mouseHOR_Sphase.txt')
g2m.genes<-readLines('/Users/carloleonardi/Desktop/R/cool-salmon/cell_cycle_genes_mouseHOR_G2Mphase.txt')
pbmc<-CellCycleScoring(pbmc, g2m.features=g2m.genes[g2m.genes %in% rownames(pbmc@assays$RNA@data)], s.features=s.genes[s.genes %in% rownames(pbmc@assays$RNA@data)], set.ident = FALSE)
pbmc@meta.data$CC.Difference <- pbmc@meta.data$S.Score - pbmc@meta.data$G2M.Score

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

pbmc <- SCTransform(pbmc, vars.to.regress = c("CC.Difference","percent.mt"), verbose = FALSE, return.only.var.genes=FALSE)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = c(0.2,0.5,1))

pbmc <- RunUMAP(pbmc, dims = 1:15)
pbmc <- RunTSNE(pbmc, dims = 1:15)

tsne_highperpl <- RunTSNE(pbmc, dims = 1:15, perplexity = 120, tsne.method="FIt-SNE")

DimPlot(tsne_highperpl, reduction = "tsne")








#mean expression di 4 o più geni (signature, scegli marker biologicamente significativi)
#umap, colori solo un cluster
#test su risultati diversi --> cfr tra log fc, definisci deg con filtri

ggplot(metadata,aes(x=data_expr,y=seurat_clusters))
+geom_density_ridges()