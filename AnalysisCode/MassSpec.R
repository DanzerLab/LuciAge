library(qs)
library(DEP)
library(limma)
library(readxl)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(pheatmap)
library(stringr)
library(ComplexHeatmap)
library(GOSemSim)
library(simplifyEnrichment)
library(rrvgo)
library(wordcloud)
library(pathfindR)
library(tidyr)
library(radiant.data)
library(ggrepel)
library(writexl)

ms_1 <- read_excel("path to MassSpec MaxLFQ raw data 1")
ms_2 <- read_excel("path to MassSpec MaxLFQ raw data 2")
metadata_1 <- data.frame(
  sample = c("3198", "3200", "3250", "3254", 
             "3206", "3207", "3211","3244",
             "2022", "2048", "2073", "2289",
             "2159","2152","2216","2257",
             "2038","2120","2453","2454",
             "1762","1764","1770","1779",
             "1756","1758","1854","2195",
             "1879","1887","1882","1999",
             "3233","3237"),
  sex = c("male","female", "male","female",
          "female","female", "male" ,"male",
          "male","male", "male" ,"female",
          "female","male", "male" ,"male",
          "female","female", "male" ,"male",
          "female","female", "female" ,"female",
          "female","female", "female" ,"male",
          "female","female", "male" ,"male",
          "male", "female"),
  age = c("6M","6M","6M","6M",
          "6M","6M","6M","6M",
          "16M","16M","16M","16M",
          "16M","16M","16M","16M",
          "16M","16M","16M","16M",
          "24M","24M","24M","24M",
          "24M","24M","24M","24M",
          "24M","24M","24M","24M",
          "6M","6M"),
  lucitrt = c("ON","ON","ON","ON",
              "OFF","OFF","OFF","OFF",
              "lateON","lateON","lateON","lateON",
              "ON","ON","ON","ON",
              "OFF","OFF","OFF","OFF",
              "lateON","lateON","lateON","lateON",
              "ON","ON","ON","ON",
              "OFF","OFF","OFF","OFF",
              "ON_Casin", "ON_Casin"),
  batch = c("Batch1", "Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1","Batch1","Batch1",
            "Batch1","Batch1")
)

metadata_2 <- data.frame(
  sample = c("3210","3245","2367","2431", "3233","3237"),
  sex = c("female", "male","female","female","male","female"),
  age = c("6M","6M","24M","24M","6M","6M"),
  lucitrt = c("OFF","OFF","OFF","OFF","ON_Casin","ON_Casin"),
  batch = c("Batch2","Batch2","Batch2","Batch2","Batch2","Batch2")
)

metadata_1$condition <- paste0(metadata_1$lucitrt, metadata_1$age)
metadata_2$condition <- paste0(metadata_2$lucitrt, metadata_2$age)

# Data preparation ----
# Create vector of columns to keep
lfq_columns_to_keep_1 <- paste("LFQ intensity", metadata_1$sample, sep=" ")
lfq_columns_to_keep_2 <- paste("LFQ intensity", metadata_2$sample, sep=" ")
non_lfq_columns_1 <- names(ms_1)[!grepl("^LFQ intensity", names(ms_1))]
non_lfq_columns_2 <- names(ms_2)[!grepl("^LFQ intensity", names(ms_2))]
columns_to_keep_1 <- c(non_lfq_columns_1, lfq_columns_to_keep_1)
columns_to_keep_2 <- c(non_lfq_columns_2, lfq_columns_to_keep_2)

# Filter the ms_data dataframe to only include these columns
data_1 <- ms_1[, columns_to_keep_1, drop = FALSE]
data_2 <- ms_2[, columns_to_keep_2, drop = FALSE]
# change column names
colnames(data_1) <- gsub("T: ", "", colnames(data_1)) # Remove "T: " and "N: "
colnames(data_1) <- gsub("N: ", "", colnames(data_1)) # Remove "T: " and "N: "
colnames(data_1) <- gsub(" \\+ ", "_", colnames(data_1)) # Replace " + " with "_"
colnames(data_1) <- gsub(" ", ".", colnames(data_1)) # Replace spaces with dots
colnames(data_2) <- gsub("T: ", "", colnames(data_2)) # Remove "T: " and "N: "
colnames(data_2) <- gsub("N: ", "", colnames(data_2)) # Remove "T: " and "N: "
colnames(data_2) <- gsub(" \\+ ", "_", colnames(data_2)) # Replace " + " with "_"
colnames(data_2) <- gsub(" ", ".", colnames(data_2)) # Replace spaces with dots

data_1$Gene.names %>% duplicated() %>% any()
data_1$Protein.IDs%>% duplicated() %>% any()
data_2$Gene.names %>% duplicated() %>% any()
data_2$Protein.IDs%>% duplicated() %>% any()
# Make a table of duplicated gene names
data_unique_1 <- make_unique(data_1, "Gene.names", "Protein.IDs", delim = ";")
data_unique_1$name %>% duplicated() %>% any()
data_unique_2 <- make_unique(data_2, "Gene.names", "Protein.IDs", delim = ";")
data_unique_2$name %>% duplicated() %>% any()

colnames(data_unique_1) <- paste0(colnames(data_unique_1),"-1")
colnames(data_unique_2) <- paste0(colnames(data_unique_2),"-2")
# creating SummarizedExperiment
LFQ_columns_1 <- grep("LFQ.", colnames(data_unique_1)) # get LFQ column numbers
LFQ_columns_2 <- grep("LFQ.", colnames(data_unique_2)) # get LFQ column numbers
lfq_data_1 <- data_unique_1[, LFQ_columns_1]
lfq_data_2 <- data_unique_2[, LFQ_columns_2]


metadata_1$sample <- paste0(metadata_1$sample, "-1")
metadata_2$sample <- paste0(metadata_2$sample, "-2")

experimental_design_1 <- metadata_1[match(gsub("LFQ.intensity.", "", colnames(lfq_data_1)), metadata_1$sample), ]
experimental_design_2 <- metadata_2[match(gsub("LFQ.intensity.", "", colnames(lfq_data_2)), metadata_2$sample), ]
colnames(experimental_design_1) <- c("label","sex","age","lucitrt","batch","condition")
colnames(experimental_design_2) <- c("label","sex","age","lucitrt","batch","condition")
experimental_design_1$sample <- sub("-\\d+$", "",experimental_design_1$label)
experimental_design_2$sample <- sub("-\\d+$", "",experimental_design_2$label)
experimental_design_1$replicate <- c(rep(c(1:4),8),c(1,2))
experimental_design_2$replicate <- c(rep(c(5:6),3))

# merge the metadata and the counts
merged_experimental_design <- rbind(experimental_design_1, experimental_design_2)

data_unique_1$name <- data_unique_1$`name-1`
data_unique_2$name <- data_unique_2$`name-2`
data_unique_1$ID <- data_unique_1$`ID-1`
data_unique_2$ID <- data_unique_2$`ID-2`

data_unique_1$`name-1` <- NULL
data_unique_2$`name-2` <- NULL
data_unique_1$`ID-1` <- NULL
data_unique_2$`ID-2` <- NULL

data_unique_1$identifier <- paste0(data_unique_1$name,"-",data_unique_1$ID)
data_unique_2$identifier <- paste0(data_unique_2$name,"-",data_unique_2$ID)

IDs_to_keep <- intersect(data_unique_1$identifier,data_unique_2$identifier) 

data_unique_1 <- data_unique_1[which(data_unique_1$identifier %in% IDs_to_keep),]
data_unique_2 <- data_unique_2[which(data_unique_2$identifier %in% IDs_to_keep),]

merged_data <- merge(data_unique_1, data_unique_2, by = c("ID","name","identifier"), all = TRUE)
LFQ_columns <- grep("LFQ.", colnames(merged_data)) # get LFQ column numbers

data_se <- make_se(merged_data, LFQ_columns, merged_experimental_design)
data_unique <- merged_data

intConditions <- c("ON6M_1","ON6M_2","ON6M_3","ON6M_4",
  "OFF6M_3","OFF6M_4","OFF6M_5","OFF6M_6",
  "ON16M_1","ON16M_2","ON16M_3","ON16M_4",
  "OFF16M_1","OFF16M_2","OFF16M_3","OFF16M_4",
  "ON24M_1","ON24M_2","ON24M_3","ON24M_4",
  "OFF24M_3","OFF24M_4","OFF24M_1","OFF24M_2",
  "lateON16M_1","lateON16M_2","lateON16M_3","lateON16M_4",
  "OFF16M_1","OFF16M_2","OFF16M_3","OFF16M_4",
  "lateON24M_1","lateON24M_2","lateON24M_3","lateON24M_4"
  )

assay_data <- as.data.frame(assay(data_se))
assay_data <- rownames_to_column(assay_data, var = "name")
assay_data$ID <- data_unique$ID
## change intConditions and run this to subset only the conditions for the comparison
colData <- as.data.frame(colData(data_se))
colData <- colData[which(paste0(colData$condition, "_", colData$replicate) %in% intConditions),]
colData$label <- colData$ID
data_se <- make_se(assay_data, which(colnames(assay_data) %in% intConditions), colData)
## starting QC
plot_frequency(data_se)
### filtering out batch detection proteins
assay_data <- assay(data_se) 
assay_data[is.na(assay_data)] <- 0
detection_matrix <- assay_data > 0
col_names <- colnames(assay_data)
batch1_cols <- grep("_1$|_2$|_3$|_4$", col_names)
batch2_cols <- grep("_5$|_6$", col_names)
valid_proteins <- (rowSums(detection_matrix[, batch1_cols]) >= 3) & (rowSums(detection_matrix[, batch2_cols]) >= 2)
data_se_filtered <- data_se[valid_proteins, ]

data_filt <- filter_missval(data_se_filtered, thr = 1) # only keeping proteins that are in 3/4 replicates of at least one condition
#data_filt <- filter_proteins(data_se, "complete") # only keeping proteins that are in all samples

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_se)+ geom_text(aes(label = sprintf("%.1f", sum), y= sum),  
                                   vjust = 3) 
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Visualize normalization
plot_normalization(data_filt, data_norm)
meanSdPlot(data_filt) 
meanSdPlot(data_norm) 

# mixed imputation according to technical or biological protein effect ----
plot_missval(data_norm)
na_summary <- get_df_long(data_norm) %>%
  group_by(name, condition) %>%
  summarize(NAs = sum(is.na(intensity)))

proteins_MNAR <- na_summary %>%
  dplyr::filter(NAs >= 3) %>%
  pull(name) %>%
  unique()

MNAR <- names(data_norm) %in% proteins_MNAR
data_imp <- impute(
  data_norm, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "MinProb") # imputation function for MNAR

plot_detect(data_norm)
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

p.pca <- plot_pca(data_imp, x = 1, y = 2, indicate = c("condition", "batch"), 
                  n = 500, # use all detected proteins
                  point_size = 4, label=F)
p.pca

mat <- assay(data_imp)
mm <- model.matrix(~condition, colData(data_imp))
mat <- limma::removeBatchEffect(mat, batch1=data_imp$sex, batch2 = data_imp$batch, design=mm)
data_imp_noBatch <- data_imp
assay(data_imp_noBatch) <- mat
plot_pca(data_imp_noBatch, x = 1, y = 2, indicate = c("condition", "batch"), 
         n = 500, # use all detected proteins
         point_size = 4, label=F)

cor_matrix <- plot_cor(data_imp_noBatch, 
                       significant = F, 
                       lower = 0, 
                       upper = 1, 
                       pal = "GnBu",
                       indicate = c("condition", "replicate"), 
                       plot = F)
pheatmap(cor_matrix)

data_diff_all_contrasts <- test_diff(data_imp, type = "manual",
                                     test = c("ON6M_vs_OFF6M","ON16M_vs_OFF16M", "lateON16M_vs_OFF16M",
                                              "ON24M_vs_OFF24M", "lateON24M_vs_OFF24M", "OFF24M_vs_OFF16M",
                                              "OFF24M_vs_OFF6M", "OFF16M_vs_OFF6M"),
                                     design_formula = formula(~0+condition+sex+batch))

dep <- add_rejections(data_diff_all_contrasts, alpha = 0.1, lfc = 0) #differential object