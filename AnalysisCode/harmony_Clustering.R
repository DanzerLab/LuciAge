library(harmony)
library(Seurat)
library(qs)

seurat_integrated <- *seurat object after all filtering has been done post scvi and clustering; see scVI_Clustering.txt*
seurat_integrated <- SCTransform(seurat_integrated, vars.to.regress = c("percent.mt","percent.ribo", "batchCorrect"), assay = "RNA",
                                 min_cells = 3, variable.features.n = 4000, method = "glmGamPoi", return.only.var.genes = TRUE)
seurat_integrated <- RunPCA(seurat_integrated)
seurat_integrated <- RunHarmony(seurat_integrated, 
                                group.by.vars = c("batchCorrect"), 
                                assay.use = "SCT", reduction.save = "harmony")

resolutions_list <- c(0.1,0.2,0.4,0.6)

seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:30, reduction = "harmony")
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = resolutions_list)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30, reduction = "harmony", n.components = 2, 
                             n.neighbors = 200, min.dist = 0.3, seed.use = 12345) # n.neighbors and min.dist can be changed. This is just for visualization