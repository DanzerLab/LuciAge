library(Matrix)
library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(qs)
library(SoupX)

# replace any *...* with respective text

# Opening files and removing background mRNA using SoupX ----
samples <- paste0("Syn", c(1:17)) # lucicompen
for (sample in samples) {
  print(paste0("Processing Sample ", sample))
  sc1 = load10X(paste0("*path to directory with all samples*", sample,"/outs/"))
  if (sample %in% c("Syn10","Syn4")){
    sc = autoEstCont(sc1, tfidfMin = 1.3) # tfidfMin is default 1. might have to increase it 
    out = adjustCounts(sc, roundToInt = TRUE)
    assign(sample, CreateSeuratObject(out))
    success = TRUE
  } else {
    success = FALSE
    tfidfMin = 1.0
    # It can happen that there is no or neglectable amount of ambient mRNA in the sample. We defined that point as 0.9 of tfidfMin
    # Meaning if a sample does not find enough "Soup Genes" at a tfidfMin of 0.9 or larger we will just work without SoupX for this sample
    while (!success && tfidfMin > 0.9) {
      tryCatch({
        print(tfidfMin)
        sc = autoEstCont(sc1, tfidfMin = tfidfMin)
        success = TRUE  # If autoEstCont succeeds, exit the loop
        out = adjustCounts(sc, roundToInt = TRUE)
        assign(sample, CreateSeuratObject(out))
      }, error = function(e) {
        tfidfMin <<- tfidfMin - 0.01  # Decrease tfidfMin by 0.01 if there's an error
      })
    }
  }
  if (!success) {
    matrix_dir <- paste0("*path to directory with all samples*", sample, "/outs/filtered_feature_bc_matrix/")
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- Read10X(data.dir = matrix_dir)
    
    assign(sample, CreateSeuratObject(counts = mat, project = paste0("Luci_mice_", sample), min.cells = 5, min.features = 100))
  }
}

merged_seurat_SoupX <- merge(x=Syn1 3, y=c(Syn2, Syn3, 
                                         Syn4, Syn5, Syn6, Syn7, Syn8, 
                                         Syn9, Syn10, Syn11, Syn12, Syn13, 
                                         Syn14, Syn15, Syn16, Syn17), 
                             add.cell.id=c("Syn1_24M_20ON","Syn2_24M_ON","Syn3_24M_OFF",
                                           "Syn4_16M_ON","Syn5_16M_ON","Syn6_16M_OFF","Syn7_24M_OFF","Syn8_16M_OFF",
                                           "Syn9_16M_12ON","Syn10_16M_12ON","Syn11_6M_ON", "Syn12_6M_OFF", "Syn13_6M_ON", 
                                           "Syn14_6M_OFF", "Syn15_24M_20ON", "Syn16_24M_ON", "Syn17_24M_ON"))
merged_seurat_SoupX <- JoinLayers(merged_seurat_SoupX)

qsave(merged_seurat_SoupX, file = "*path for saving*", nthreads = 4)
