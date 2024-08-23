# classify Sex using cellXY
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(cellXY)
library(CellBench)
library(BiocStyle)
library(scater)
library(qs)
library(Seurat)
library(rlang) # might have to reinstall
library(caret) # might have to reinstall
# only get samples that have male,female
filtered_seurat <- qread(file = "*path to filtered seurat; see QualityControl.R*")
filtered_seurat$sampleDeMulti <- filtered_seurat$sample
rel_samples <- unique(filtered_seurat$sample[which(filtered_seurat$sex == "male,female")])
for (s in rel_samples){
  seurat_sex <- subset(filtered_seurat, sample == s)
  DefaultAssay(seurat_sex) <- "RNA"
  counts_matrix <- as.matrix(GetAssayData(seurat_sex, "RNA"))
  metadata <- seurat_sex@meta.data
  ## sce; make sure rows are genes and columns are cells, and check if metadata is correct
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix), colData = metadata)
  counts <- counts(sce)
  ## finding sex doublets
  doublets <- findMfDoublet(counts, genome = "Mm")
  
  singlets <- counts[,doublets$prediction=="Singlet"]
  sex <- classifySex(singlets, genome="Mm",qc=TRUE)
  
  metadata <- filtered_seurat@meta.data
  male <- rownames(subset(sex, prediction == "Male"))
  female <- rownames(subset(sex, prediction == "Female"))
  ### changing male,female sample to either female or male
  metadata$sex[which(metadata$cells %in% male)] <- "male"
  metadata$sex[which(metadata$cells %in% female)] <- "female"
  metadata$sampleDeMulti[which(metadata$cells %in% male)] <- paste0(seurat_sex$sample[1],"_male")
  metadata$sampleDeMulti[which(metadata$cells %in% female)] <- paste0(seurat_sex$sample[1],"_female")
  filtered_seurat@meta.data <- metadata
}

filtered_seurat <- subset(filtered_seurat, sex == "male,female", invert = TRUE)
metadata <- filtered_seurat@meta.data
metadata$sampleDeMulti[which(metadata$sex == "male")] <- paste0(metadata$sample[which(metadata$sex == "male")], "_male")
metadata$sampleDeMulti[which(metadata$sex == "female")] <- paste0(metadata$sample[which(metadata$sex == "female")], "_female")
filtered_seurat@meta.data <- metadata
qsave(filtered_seurat, file = "*path to save*", nthreads = 4)
