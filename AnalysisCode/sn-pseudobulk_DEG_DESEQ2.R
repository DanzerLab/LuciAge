library(qs)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(stringr)
library(cowplot)
library(DESeq2)
library(scran)
library(pheatmap)
library(ggrepel)
library(radiant.data)
library(writexl)
library(readxl)
library(BiocParallel)
library(AnnotationHub)
library(fdrtool)
library(locfdr)
library(magrittr)
library(sva)
library(hms)

# replace *...* with respective text

# final code for for loop with sva, no ribo, no mito, no sex, balanced 1.5 ----
seurat_integrated <- qread(file = "*final filtered, integrated and annotated seurat object*")

ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"))
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

comps <- c("synON16M_vs_synOFF16M","synlateON16M_vs_synOFF16M", "synON24M_vs_synOFF24M", "synlateON24M_vs_synOFF24M", "synON6M_vs_synOFF6M", 
           "synOFF16M_vs_synOFF6M", "synOFF24M_vs_synOFF6M", "synOFF24M_vs_synOFF16M")
comparison_to_conditions <- list(
  "synON16M_vs_synOFF16M" = c("synON16M", "synOFF16M"),
  "synlateON16M_vs_synOFF16M" = c("synlateON16M", "synOFF16M"),
  "synON24M_vs_synOFF24M" = c("synON24M", "synOFF24M"),
  "synlateON24M_vs_synOFF24M" = c("synlateON24M", "synOFF24M"),
  "synON6M_vs_synOFF6M" = c("synON6M", "synOFF6M"),
  "synOFF16M_vs_synOFF6M" = c("synOFF16M", "synOFF6M"), 
  "synOFF24M_vs_synOFF6M"= c("synOFF24M", "synOFF6M"), 
  "synOFF24M_vs_synOFF16M" = c("synOFF24M", "synOFF16M")
)

number_of_loops <- 1
Idents(seurat_integrated) <- "CellLevel3"

cellN <- as.data.frame(table(seurat_integrated@meta.data$condition))
lucitrt <- unique(seurat_integrated@meta.data$condition)


seurat_integrated$bulk <- paste0(seurat_integrated$sample,"-",seurat_integrated$condition,"-", seurat_integrated$age,"-", 
                                 "-",seurat_integrated$batch, "-", seurat_integrated$sex)

ctName_vector <- as.character(unique(Idents(seurat_integrated)))

df_results_sva <- data.frame(
  comparison = character(),
  celltype = character(),
  ON24 = integer(),
  ON16 = integer(),
  OFF24 = integer(),
  OFF16 = integer(),
  ON20 = integer(),
  ON12 = integer(),
  ON6 = integer(),
  OFF6 = integer(),
  nUMI = integer(),
  Nsigres = integer(),
  sigres_perct = numeric(),
  downregulated = integer(),
  upregulated = integer(),
  downregulated_pct = numeric(),
  upregulated_pct = numeric(),
  stringsAsFactors = FALSE
)
sigres_list_sva <- list()
allres_list_sva <- list()
progress = 0
start_time <- Sys.time()
error_messages <- character()
for (ctName in ctName_vector){
  tryCatch({
    print(paste0("Differential for: ", ctName))
    seurat_deg <- subset(seurat_integrated, idents = ctName)
    for (i  in c(1:number_of_loops)){
      progress <- progress + 1
      print(paste0("Progress: ", progress,"/",number_of_loops*length(ctName_vector)*length(comps)))
      
      current_time <- Sys.time()
      duration <- as_hms(difftime(current_time, start_time))
      print(paste0("Code has been running for ", duration, "  HH:MM:SS"))
      # add print for time statement. How long it was running
      
      cellN <- as.data.frame(table(seurat_deg@meta.data$condition))
      for (comp in comps){
        print(paste0("Processing: ",comp))
        conditions <- comparison_to_conditions[[comp]]
        comp_index <- which(comps == comp)
        ## checking for batch effect withing conditions
        ###Condition 1
        if (length(unique(seurat_deg$batch[which(seurat_deg$condition %in% conditions[1])])) > 1){
          print(paste0("Subsampling ", conditions[1]))
          seurat_deg_cond <- subset(seurat_deg, condition %in% conditions[1])
          number_of_cells <- round(min(table(seurat_deg_cond$sample)) * 1.5)
          keep_cells_list <- list()
          ### Go through every sample and check if some have to be subsampled
          for(s in unique(seurat_deg_cond$sample)){
            if(length(seurat_deg_cond$sample[which(seurat_deg_cond$sample == s)]) > number_of_cells){
              sample_seurat <- subset(seurat_deg_cond, sample %in% s)
              set.seed(23193)
              sub_cells <- sample(sample_seurat$cells,number_of_cells)
              keep_cells_list[[s]] <- sub_cells
            }else{
              keep_cells_list[[s]] <- colnames(seurat_deg_cond[,which(seurat_deg_cond$sample == s)])
            }
          }
          batch_keep_cells_cond1 <- c(unlist(keep_cells_list, use.names = FALSE))
          print(table(seurat_deg_cond$sample))
          print(number_of_cells)
          print(table(seurat_deg_cond$sample[which(seurat_deg_cond$cells %in% batch_keep_cells_cond1)]))
        }else{
          batch_keep_cells_cond1 <- seurat_deg$cells[which(seurat_deg$condition == conditions[1])]
        }
        ## condition 2
        if (length(unique(seurat_deg$batch[which(seurat_deg$condition %in% conditions[2])])) > 1){
          print(paste0("Subsampling ", conditions[2]))
          seurat_deg_cond <- subset(seurat_deg, condition %in% conditions[2])
          number_of_cells <- round(min(table(seurat_deg_cond$sample)) * 1.5)
          keep_cells_list <- list()
          ### Go through every sample and check if some have to be subsampled
          for(s in unique(seurat_deg_cond$sample)){
            if(length(seurat_deg_cond$sample[which(seurat_deg_cond$sample == s)]) > number_of_cells){
              sample_seurat <- subset(seurat_deg_cond, sample %in% s)
              set.seed(23193)
              sub_cells <- sample(sample_seurat$cells,number_of_cells)
              keep_cells_list[[s]] <- sub_cells
            }else{
              keep_cells_list[[s]] <- colnames(seurat_deg_cond[,which(seurat_deg_cond$sample == s)])
            }
          }
          batch_keep_cells_cond2 <- c(unlist(keep_cells_list, use.names = FALSE))
          print(table(seurat_deg_cond$sample))
          print(number_of_cells)
          print(table(seurat_deg_cond$sample[which(seurat_deg_cond$cells %in% batch_keep_cells_cond2)]))
        }else{
          batch_keep_cells_cond2 <- seurat_deg$cells[which(seurat_deg$condition == conditions[2])]
        }
        
        ## now balancing between comparisons
        min_final <- round(min(length(batch_keep_cells_cond1), length(batch_keep_cells_cond2)) * 1.5)
        if (comp == "synON6M_vs_synOFF6M"){
          min_final <- 326
        }
        
        if (length(batch_keep_cells_cond1) > min_final){
          set.seed(23193)
          batch_keep_cells_cond1_final <- sample(batch_keep_cells_cond1, min_final)
        }else{
          batch_keep_cells_cond1_final <- batch_keep_cells_cond1
        }
        
        if (length(batch_keep_cells_cond2) > min_final){
          set.seed(23193)
          batch_keep_cells_cond2_final <- sample(batch_keep_cells_cond2, min_final)
        }else{
          batch_keep_cells_cond2_final <- batch_keep_cells_cond2
        }
        cellN$Freq[cellN$Var1 == conditions[1]] <- as.numeric(length(batch_keep_cells_cond1_final))
        cellN$Freq[cellN$Var1 == conditions[2]] <- as.numeric(length(batch_keep_cells_cond2_final))
        
        final_keep_cells <- c(batch_keep_cells_cond1_final, batch_keep_cells_cond2_final)
        seurat_subset <- subset(seurat_deg, cells %in% final_keep_cells)
        print("Final Cellnumbers: ")
        print(table(seurat_subset$sample))
        print(table(seurat_subset$condition))
        # Extract counts
        counts <- AggregateExpression(
          object = seurat_subset,
          group.by = "bulk",  # or whatever grouping variable you have
          assay = "RNA",
          slot = "counts",
          return.seurat = FALSE
        )
        counts <- counts$RNA
        keep <- rowSums(counts>3) >= max(length(unique(seurat_deg$sample[which(seurat_deg$condition == conditions[1])])),length(unique(seurat_deg$sample[which(seurat_deg$condition == conditions[2])])))
        counts <- counts[keep,]
        counts <- counts[!grepl("^mt-", rownames(counts)), ]
        counts <- counts[-which(rownames(counts) %in% rb), ]
        counts <- counts[-which(rownames(counts) %in% c("Tsix","Xist")), ]
        ## to minimize potential contamination of results from batch effect we removed genes that were differentially expressed in batch4 vs batch3 for comparisons
        ## that have confounding design with batch4 vs batch3
        ## This is a very conservative approach minimizing false positive results by risking loss of some true positives
        if (comp %in% c("synOFF24M_vs_synOFF6M", "synOFF16M_vs_synOFF6M")){
          counts <- counts[-which(rownames(counts) %in% c("Camk1d","Cdk8","Lars2","Gm42418","AY036118","Tshz2","Gphn","Ptgds","Cacnb2","Pdzrn3","Auts2",
                                                          "2900097C17Rik", "Kirrel3","Pde4d","Gm13052","Adarb2","Nav3","4930555F03Rik",
                                                          "Nrg1","Hnrnpa1","Tox","Nrg3","Nicn1","Gm49678","Rgs9","Mast4","Eif2s3y","Cdh20",
                                                          "4930445B16Rik", "Frmd4b","Kcnip4")), ]
        }
        ## edit metadata for bulk
        colData <- data.frame(bulk = colnames(counts))
        
        rownames(colData) <- colnames(counts)
        
        colData$age <- NA
        colData$age[which(str_detect(colData$bulk, ".6M"))] <- "6M"
        colData$age[which(str_detect(colData$bulk, ".16M"))] <- "16M"
        colData$age[which(str_detect(colData$bulk, ".24M"))] <- "24M"
        
        colData$sex <- NA
        colData$sex[which(str_detect(colData$bulk, ".male"))] <- "male"
        colData$sex[which(str_detect(colData$bulk, ".female"))] <- "female"
        colData$sex[which(str_detect(colData$bulk, ".female.male"))] <- "female.male"
        
        colData$batch <- NA
        colData$batch[which(str_detect(colData$bulk, ".Batch3"))] <- "Batch3"
        colData$batch[which(str_detect(colData$bulk, ".Batch4"))] <- "Batch4"
        
        colData$lucitrt <- "NA"
        colData$lucitrt[which(str_detect(colData$bulk, ".synON24M."))] <- "synON24M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synON6M."))] <- "synON6M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synON16M."))] <- "synON16M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synOFF6M."))] <- "synOFF6M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synOFF24M."))] <- "synOFF24M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synOFF16M."))] <- "synOFF16M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synlateON24M."))] <- "synlateON24M"
        colData$lucitrt[which(str_detect(colData$bulk, ".synlateON16M."))] <- "synlateON16M"
        
        colData[] <- lapply(colData, as.factor)
        
        dds <- DESeqDataSetFromMatrix(counts, 
                                      colData = colData, 
                                      design =  ~lucitrt)
        reference_level <- conditions[2]
        dds$lucitrt <- factor(colData$lucitrt, levels = c(reference_level, setdiff(unique(colData$lucitrt), reference_level)))
        
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        
        norm.cts <- counts(dds, normalized=TRUE)
        mm <- model.matrix(~ lucitrt, colData(dds))
        mm0 <- model.matrix(~ 1, colData(dds))
        norm.cts <- norm.cts[rowSums(norm.cts) > 10,]
        fit <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv = 1)
        dds$SV1 <- fit$sv[,1]
        design(dds) <- ~ SV1 + lucitrt 
        tryCatch({
          fit <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv = 2)
          dds$SV1 <- fit$sv[,1]
          dds$SV2 <- fit$sv[,2]
          design(dds) <- ~ SV1 + SV2 + lucitrt}, 
          error = function(e){
            print("Using one SV")
          })
        
        dds <- DESeq(dds, test = "Wald")
        contrast <- resultsNames(dds)[length(resultsNames(dds))]
        res <- results(dds, alpha = 0.05, name = contrast)
        FDR.ddsRes <- fdrtool(res$stat, statistic= "normal", plot = T)
        res2 <- lfcShrink(dds, res = res, type = "normal", coef = contrast)
        res2$pvalue  <- FDR.ddsRes$pval
        res2$padj <- p.adjust(res2$pvalue, method = "BH")
        res_tbl_bulk <-as.data.frame(res2)
        res_tbl_bulk$gene <- rownames(res_tbl_bulk)
        res_tbl_bulk$stat <- res$stat
        sig_res <- dplyr::filter(res_tbl_bulk, padj < 0.05) %>% dplyr::arrange(padj)
        if (ctName %in% c("GABA-Rmst/Tshz2 Neuron", "GABA-Ano1/2 Neuron")){
          ctName <- gsub("/", "_", ctName)
        }
        sigres_name <- paste0(ctName,"-",comp,"-",i)
        
        sigres_list_sva[[sigres_name]] <- sig_res
        allres_list_sva[[sigres_name]] <- res_tbl_bulk
        
        df_row <- data.frame(
          comparison = paste0(comp,"-",i),
          celltype = ctName,
          ON24 = ifelse("synON24M" %in% conditions, cellN[cellN$Var1 == "synON24M", "Freq"], 0),
          ON16 = ifelse("synON16M" %in% conditions, cellN[cellN$Var1 == "synON16M", "Freq"], 0),
          OFF24 = ifelse("synOFF24M" %in% conditions, cellN[cellN$Var1 == "synOFF24M", "Freq"], 0),
          OFF16 = ifelse("synOFF16M" %in% conditions, cellN[cellN$Var1 == "synOFF16M", "Freq"], 0),
          ON6 = ifelse("synON6M" %in% conditions, cellN[cellN$Var1 == "synON6M", "Freq"], 0),
          OFF6 = ifelse("synOFF6M" %in% conditions, cellN[cellN$Var1 == "synOFF6M", "Freq"], 0),
          ON20 = ifelse("synlateON24M" %in% conditions, cellN[cellN$Var1 == "synlateON24M", "Freq"], 0),
          ON12 = ifelse("synlateON16M" %in% conditions, cellN[cellN$Var1 == "synlateON16M", "Freq"], 0),
          nUMI = sum(colSums(counts(dds))),
          Nsigres = nrow(sig_res),
          sigres_perct = length(sig_res$gene) / length(res_tbl_bulk$gene),
          downregulated = sum(sig_res$log2FoldChange < 0),
          upregulated = sum(sig_res$log2FoldChange > 0),
          downregulated_pct = (sum(sig_res$log2FoldChange < 0) / nrow(sig_res)) * 100,
          upregulated_pct = (sum(sig_res$log2FoldChange > 0) / nrow(sig_res)) * 100,
          stringsAsFactors = FALSE
        )
        df_results_sva <- rbind(df_results_sva, df_row)
        print(paste0("Added Information Data Frame: "))
        print(head(df_row))
        
        write_xlsx(sig_res, paste0("*path to save results*",ctName,"_",comp,"_sigres.xlsx"))
        write_xlsx(res_tbl_bulk, paste0("*path to save results*",ctName, "_",comp,"_allres.xlsx"))
      }
    }
  },error = function(e) {
    # Save the error message with the comparison identifier
    error_message <- sprintf("Error in %s: %s", condition, e$message)
    print(paste0("error_message: ", error_message))
    message(error_message)
    # Append the error message to the vector
    error_messages <<- c(error_messages,error_message)
  })
}
write_xlsx(df_results_sva,"*path to save results*/infoDF_all1.5.xlsx")