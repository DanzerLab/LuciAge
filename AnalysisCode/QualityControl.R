library(qs)
library(scDblFinder)
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(stringr)
library(AnnotationHub)
library(ensembldb)
library(gridExtra)
library(RCurl)
library(ggplot2)
library(cowplot)

# Replace all *...* with respective text

# open merged_seurat file. Check path ----
merged_seurat <- qread(file = "*path of mergedSeurat; see SoupX_MergedSeurat.R*")

ah <- AnnotationHub()
ahDb <- AnnotationHub::query(ah, pattern = c("Mus musculus", "EnsDb"))
#getting most recent ID
id <- ahDb %>%  
  mcols() %>%
  rownames() %>%
  tail(n = 1)
edb <- ah[[id]] #download ensembl database
annotations <- genes(edb, return.type = "data.frame") # Extract gene-level information from database
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
# Extract IDs for mitochondrial genes
mt <- annotations %>%
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_name)
# Extract IDs for ribosomal genes
rb <- annotations %>%
  dplyr::filter(grepl("ribosomal", description, ignore.case = TRUE)) %>%
  dplyr::pull(gene_name)
# calculating percent.mt and percent.ribo
counts <- GetAssayData(merged_seurat, slot = "counts")
mtUMI <- Matrix::colSums(counts[which(rownames(counts) %in% mt),], na.rm = T)
rbUMI <- Matrix::colSums(counts[which(rownames(counts) %in% rb),], na.rm = T)
mtperUMI <- mtUMI/merged_seurat$nCount_RNA
rbperUMI <- rbUMI/merged_seurat$nCount_RNA
merged_seurat[["percent.mt"]] <- mtperUMI
merged_seurat[["percent.ribo"]] <- rbperUMI
merged_seurat[["log10GenesPerUMI"]] <- log10(merged_seurat$nFeature_RNA)/log10(merged_seurat$nCount_RNA)
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
## add extra information into metadata
metadata$sample <- NA
for (i in c(1:17)) {
  sample_name <- paste0("Syn", i)
  metadata$sample[which(str_detect(metadata$cells, paste0("^", sample_name, "_")))] <- sample_name
}

metadata$age <- NA
metadata$age[which(str_detect(metadata$cells, "_6M_"))] <- "6M"
metadata$age[which(str_detect(metadata$cells, "_24M_"))] <- "24M"
metadata$age[which(str_detect(metadata$cells, "_16M_"))] <- "16M"


metadata$sex <- NA
for (i in c(1,3,10, 15,16,17)) {
  sample_name <- paste0("Syn", i)
  metadata$sex[which(str_detect(metadata$cells, paste0("^", sample_name, "_")))] <- "female"
}
for (i in c(2,4,6:8,11:14)) {
  sample_name <- paste0("Syn", i)
  metadata$sex[which(str_detect(metadata$cells, paste0("^", sample_name, "_")))] <- "female.male"
}
for (i in c(5,9)) {
  sample_name <- paste0("Syn", i)
  metadata$sex[which(str_detect(metadata$cells, paste0("^", sample_name, "_")))] <- "male"
}

metadata$nAnimals <- 2
metadata$nAnimals[which(str_detect(metadata$cells, "^Syn9_"))] <- 1
metadata$nAnimals[which(str_detect(metadata$cells, "^Syn10_"))] <- 1
metadata$nAnimals[which(str_detect(metadata$cells, "^Syn15_"))] <- 1
metadata$nAnimals[which(str_detect(metadata$cells, "^Syn16_"))] <- 1
metadata$nAnimals[which(str_detect(metadata$cells, "^Syn17_"))] <- 1

metadata$batch <- "Batch3"
metadata$batch[which(str_detect(metadata$cells, "^Syn11_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn12_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn13_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn14_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn15_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn16_"))] <- "Batch4"
metadata$batch[which(str_detect(metadata$cells, "^Syn17_"))] <- "Batch4"

metadata$lucitrt <- NA
metadata$lucitrt[which(str_detect(metadata$cells, "_OFF"))] <- "synOFF"
metadata$lucitrt[which(str_detect(metadata$cells, "_ON_"))] <- "synON"
metadata$lucitrt[which(str_detect(metadata$cells, "_20ON_"))] <- "synlateON"
metadata$lucitrt[which(str_detect(metadata$cells, "_12ON_"))] <- "synlateON"

metadata$condition <- paste0(metadata$lucitrt, metadata$age)

## add metadata back to merged_seurat
merged_seurat@meta.data <- metadata
metadata <- merged_seurat@meta.data

# filtering dataset (check the values) ----
filtered_seurat <- subset(x = merged_seurat, 
                                 subset= (nUMI >= 2500) & 
                                   (nGene >= 1500) & 
                                   (log10GenesPerUMI > 0.85) & 
                                   (percent.mt < 0.03) & (percent.ribo < 0.015))
# doublet detection using scdbl ----
## prep counts matrix and metadata for singlecellexperiment object
counts_matrix <- as.matrix(GetAssayData(filtered_seurat, "RNA"))
metadata <- filtered_seurat@meta.data
## sce; make sure rows are genes and columns are cells, and check if metadata is correct
sce <- SingleCellExperiment(assays = list(counts = counts_matrix), colData = metadata)
sce <- scDblFinder(sce, samples = "sample") ## specify samples!!!
## back to seurat. also here check if score and class assignment is correct
filtered_seurat$scDblFinder.score <- sce$scDblFinder.score
filtered_seurat$scDblFinder.class <- sce$scDblFinder.class
## subsetting singlets
filtered_seurat <- subset(x = filtered_seurat, 
                          subset= (scDblFinder.class == "singlet"))

# Save and Load seurat file ----
qsave(filtered_seurat, file = "*path to save*", nthreads = 4)

