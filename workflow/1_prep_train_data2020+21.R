## This script is for replicating the results of the 2nd winner of CMI-PB Challenge 2 (Nicky) ##
# Github repo - https://github.com/nixstix/CMI-PB2/tree/main
# Datafiles downloaded from - https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/legacy/2nd_challenge/2024-02-02/2nd_challenge/processed_datasets/training_dataset/

# REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised by the CMI-PB team. 
# 
# We will take a quick look to see what the data looks like.

## LOAD DATA

source('CMI-PB2-main/scripts/libs.R')
load("CMI-PB2-main/processed_datasets/training_dataset/train_2020+21_TRAIN.RData")

dat = readRDS('master_processed_training_data.RDS')
names(dat) # RL: The rds object has 5 dataframes - subject_specimen, abtiter_wide, plasma_cytokine_concentrations, pbmc_cell_frequency, pbmc_gene_expression


## SELECT DS FROM 2021 AS TRAINING SET
# RL: The as.character() function converts a numeric object to a string data type or a character object. If the collection is passed to it as an object, it converts all the elements of the collection to a character or string type.
dat$subject_specimen$specimen_id = as.character(dat$subject_specimen$specimen_id)
# RL: table() function is used to create a categorical representation of data with variable name and the frequency in the form of a table.
table(dat$subject_specimen$dataset)
dim(dat$subject_specimen)


# MERGE DATA
## merge all matrices to replace empty values with NAs

meta = dat$subject_specimen # RL: Extracting just the subject_specimen df from master RDS
colnames(meta) = paste("Meta.", colnames(meta), sep='') # RL: Adding the string 'Meta.' to all colnames

# ab
# RL: abtiter df in master RDS has 4 sub dfs - metadata, raw_data, normalized_data, batchCorrected_data
names(dat$abtiter_wide)
dim(dat$abtiter_wide$normalized_data)
# RL: t() is used to transpose a matrix/dataframe. as.data.frame() is used to convert an object to data frame. These objects can be Vectors, Lists, Matrices, and Factors.
x = as.data.frame(t(dat$abtiter_wide$normalized_data))
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# cytokines
# RL: plasma_cytokine_concentrations df in master RDS has 4 sub dfs - metadata, raw_data, normalized_data, batchCorrected_data
x = as.data.frame(t(dat$plasma_cytokine_concentrations$normalized_data))
colnames(x) = paste("Cytokine.", colnames(x), sep='') # RL: Adding the string 'Cytokine.' to all colnames
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))


# cell_freq
# RL: pbmc_cell_frequency df in master RDS has 4 sub dfs - metadata, raw_data, normalized_data, batchCorrected_data
x = as.data.frame(t(dat$pbmc_cell_frequency$normalized_data))
colnames(x)
colnames(x) = paste("Freq.", colnames(x), sep='') # RL: Adding the string 'Freq.' to all colnames
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# expression
# RL: pbmc_gene_expression df in master RDS has 3 sub dfs - metadata, raw_data, batchCorrected_data
## NOTE: in a normal workflow, the input for this section should not be batch-corrected or TPM data, but raw count data. If raw count data is not available, the batch-corrected data can be directly fed into MOFA without any kind of filtering or transformation.
# RL: DESeq2 expects raw counts for proper normalization, and using batch-corrected data can interfere with this process. In this script, batch-corrected data is being fed into DESeq2, which may not be ideal because DESeq2's normalization functions (like size factor estimation) are meant for raw counts, and the data has already been batch-corrected. This could distort the analysis.Ideally, you should use raw count data if available for this pipeline.
#counts = dat$pbmc_gene_expression$batchCorrected_data
counts_raw = dat$pbmc_gene_expression$raw_data
# RL: The , inside the square brackets indicates that the operation is affecting the columns (right of the comma), while all rows are included (left of the comma is empty). The which() function returns the indices/column positions where colnames(counts) == meta$Meta.specimed_id
#counts = counts[, which(colnames(counts) %in% meta$Meta.specimen_id)]
counts_raw = counts_raw[, which(colnames(counts_raw) %in% meta$Meta.specimen_id)] # RL: meta$Meta.specimen_id has 729 samples while counts_raw has 360 samples. This line filters the columns of the counts matrix based on whether the column names (samples) match the Meta.specimen_id column in the metadata. It ensures that the data used for further analysis corresponds to the samples listed in the metadata.
#counts = apply(counts,2, as.integer)
counts_raw = apply(counts_raw,2, as.integer)
# RL: Converting the column (2) datatype to integer. Syntax: apply( x, margin, function ) x: determines the input array including matrix. margin: If the margin is 1 function is applied across row, if the margin is 2 it is applied across the column. function: determines the function that is to be applied on input data.
#rownames(counts) = rownames(dat$pbmc_gene_expression$batchCorrected_data) 
rownames(counts_raw) = rownames(dat$pbmc_gene_expression$raw_data)
# setting the row names of the counts object to be the same as the row names of dat$pbmc_gene_expression$batchCorrected_data for consistency.
#counts = apply(counts,2, as.integer)
#dim(counts)
dim(counts_raw)
#write.csv(counts, file= "counts.csv", row.names = F)

# RL: This creates a DESeq2 dataset (dds) from the counts matrix and the corresponding sample metadata (colData), with a design formula specifying ~ infancy_vac. The design formula indicates that differential expression will be analyzed based on the variable infancy_vac (likely representing a grouping factor such as vaccine treatment during infancy). The line dat$pbmc_gene_expression$metadata[colnames(counts_raw), ] is selecting rows from the metadata data frame that correspond to the sample names present in counts_raw. This ensures that the metadata and the raw count data are aligned by sample, which is essential for downstream analyses like differential expression, normalization, and visualization.
table(dat$pbmc_gene_expression$metadata$infancy_vac) # Has labels aP/wP
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = dat$pbmc_gene_expression$metadata[colnames(counts) , ],
#                               design = ~ infancy_vac)
dds_from_raw <- DESeqDataSetFromMatrix(countData = counts_raw,
                              colData = dat$pbmc_gene_expression$metadata[colnames(counts_raw) , ],
                              design = ~ infancy_vac)
#dds
dds_from_raw
#dds = estimateSizeFactors(dds)
# RL: This step estimates size factors for normalization using the "median ratio method". It adjusts for differences in sequencing depth between samples.
dds_from_raw = estimateSizeFactors(dds_from_raw)
#vsd <- vst(dds, blind=TRUE)
# RL: This applies variance stabilizing transformation (VST) to the data. The blind=TRUE option means that the transformation is performed without considering the experimental design. The VST helps make the count data more comparable across genes and samples, particularly for PCA or clustering. Why VST?: VST reduces the dependence of the variance on the mean, which is important for downstream analyses like PCA. Value?: a DESeqTranform object or a matrix of transformed, normalized counts
vsd_raw <- vst(dds_from_raw, blind=TRUE)
#dim(assay(vsd))
dim(assay(vsd_raw))

#vsd$subject_id = as.factor(vsd$subject_id)
# RL: as.factor() is an R function that converts a character or integer variable into a factor. A factor in R is a data type used to represent categorical data with a fixed set of levels (i.e., distinct categories). By converting subject_id into a factor, you're telling R to treat it as categorical data rather than numeric or character data. This is particularly useful for downstream analyses, such as differential expression or visualization, where you want to group samples by subjects, categories, or conditions.
vsd_raw$subject_id = as.factor(vsd_raw$subject_id)

# RL: The terms HVF1 and HVF2 in your code likely refer to "High-Variance Factors" or "Highly Variable Factors" 1 and 2. In this context, they are used to represent the principal components (PCs) derived from Principal Component Analysis (PCA). HVF1 (or PC1): This is the first principal component, which captures the greatest variance in the data. HVF2 (or PC2): This is the second principal component, which captures the second-largest variance in the data after PC1.
# RL: These lines perform PCA (Principal Component Analysis) on the variance-stabilized data (vsd) to check how the samples cluster before removing batch effects. "dataset" and "subject_id" refer to variables in the metadata that are being used as grouping factors in the PCA plots.
# bef_noHVF1 = plotPCA(vsd, "dataset") # PCA before
bef_noHVF1_raw = plotPCA(vsd_raw, "dataset") # PCA before
#bef_noHVF2 = plotPCA(vsd, "subject_id") # PCA before
bef_noHVF2_raw = plotPCA(vsd_raw, "subject_id")

#assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))
# RL: This removes batch effects from the variance-stabilized data using the removeBatchEffect function from the limma package. The batch argument specifies the variable (in this case, dds$dataset) representing batch information. This is an additional batch correction step that helps ensure the data is not biased by batch effects when running downstream analyses like PCA.
assay(vsd_raw) <- limma::removeBatchEffect(assay(vsd_raw), batch=as.vector(dds_from_raw$dataset))
# RL: We now plot PCA after removing batch-effects to for visualizing data.
#aft_noHVF1 = plotPCA(vsd, "dataset")
aft_noHVF1_raw = plotPCA(vsd_raw, "dataset")
#aft_noHVF2 = plotPCA(vsd, "subject_id")
aft_noHVF2_raw = plotPCA(vsd_raw, "subject_id")

# mofa is being dominated by factor 1, this may be a result of poor normalisation, therefore let's limit the geneset even further, based on HVG
# RL: This identifies the top 5000 most variable genes by calculating the row variance of the variance-stabilized data and selecting the genes with the highest variance. Highly variable genes are often more informative in exploratory analyses like PCA or clustering.
#topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),5000)
topVarGenes_from_raw <- head(order(rowVars(assay(vsd_raw)),decreasing=TRUE),5000)
#topVarGenes <- rownames(vsd)[topVarGenes]
topVarGenes_from_raw <- rownames(vsd_raw)[topVarGenes_from_raw]

# RL: This filters the DESeq2 dataset to retain only the top 5000 highly variable genes. Reducing the dataset to a more informative subset can improve downstream analyses by removing noise from low-variance genes.
#dds <- dds[topVarGenes,]
dds_from_raw <- dds_from_raw[topVarGenes_from_raw,]
#dds
dds_from_raw
#dds = estimateSizeFactors(dds)
# RL: The size factors are re-estimated, and the variance stabilizing transformation is applied again, but this time only to the highly variable genes. The dimensions of the new dataset are checked.
dds_from_raw = estimateSizeFactors(dds_from_raw)
#vsd <- vst(dds, blind=TRUE)
vsd_raw <- vst(dds_from_raw, blind=TRUE)
#dim(assay(vsd))
dim(assay(vsd_raw))
# RL: This repeats the PCA analysis before and after removing batch effects but now using the highly variable gene subset. The resulting PCA plots should be more informative, focusing on the most relevant genes.
#bef_plusHVF=plotPCA(vsd, "dataset") # PCA before
bef_plusHVF_raw=plotPCA(vsd_raw, "dataset")

#assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))
assay(vsd_raw) <- limma::removeBatchEffect(assay(vsd_raw), batch=as.vector(dds_from_raw$dataset))

#aft_plusHVF=plotPCA(vsd, "dataset")
aft_plusHVF_raw=plotPCA(vsd_raw, "dataset")

#x = as.data.frame(t(assay(vsd)))
# RL: # RL: t() is used to transpose a matrix/dataframe. as.data.frame() is used to convert an object to data frame. These objects can be Vectors, Lists, Matrices, and Factors.
x = as.data.frame(t(assay(vsd_raw)))
colnames(x) = paste("Expr.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

#train_data = t(meta)
# RL: The code creates a list object train_data with different types of data (e.g., gene expression, antibody levels, cytokine levels, etc.). The data comes from the meta data frame, and each part of the list represents a different type of biological data extracted from the columns of meta. grep() finds all column names in the meta data frame that contain the particular string. t(): Transposes the data. 
train_data = list(meta = meta[,grep('Meta', colnames(meta))],
                  geneExp = t(meta[,grep('Expr.', colnames(meta))]) ,
                  ab = t(meta[,grep('^IgG', colnames(meta))]) ,
                  cytokine = t(meta[,grep('Cytokine', colnames(meta))]), 
                  cell_freq = t(meta[,grep('Freq', colnames(meta))])
)
# RL: This ensures that the columns of the various data matrices (geneExp, ab, cytokine, cell_freq) are labeled with the specimen IDs from the metadata (Meta.specimen_id). 
colnames(train_data$geneExp) = train_data$meta$Meta.specimen_id
colnames(train_data$ab) = train_data$meta$Meta.specimen_id
colnames(train_data$cytokine) = train_data$meta$Meta.specimen_id
colnames(train_data$cell_freq) = train_data$meta$Meta.specimen_id

train_data$ab[1:5, 1:13] # RL: This displays the first 5 rows and 13 columns of the antibody data to preview its structure.

lapply(train_data, dim) # RL: This applies the dim function to each element of the train_data list to check the dimensions (number of rows and columns) of each dataset.

#save(train_data, file = 'processed_datasets/training_datasets/train_2020+21.RData')

# split into test and train
# RL: The code proceeds to split the data into training and test sets based on subject IDs, ensuring that the split is done by subject rather than by individual specimens to avoid data leakage. 
idx = unique(train_data$meta$Meta.subject_id)
idx
length(idx) # RL: There are 96 subjects
set.seed(42) # RL: Sets the seed for reproducibility of the random sampling.
test.idx = sample(train_data$meta$Meta.subject_id, 0.40*length(idx),replace=FALSE) # RL: Select 40% of subjects for test set
test.idx = as.vector(train_data$meta[train_data$meta$Meta.subject_id %in% test.idx, 'Meta.specimen_id'])[[1]] # RL: Convert to specimen IDs for test set
train.idx = as.vector(train_data$meta[!train_data$meta$Meta.specimen_id %in% test.idx, 'Meta.specimen_id'])[[1]] # RL: Get remaining specimens for training set
                                                                              # RL: This constructs the training dataset (train_data2) by subsetting the original train_data list for the specimens that belong to the training set (train.idx). This list contains metadata, gene expression data, antibody data, cytokine data, and cell frequency data, all specific to the training set. train_data2: A list of the same structure as train_data, but containing only the training samples.                              
train_data2 = list(
  meta = train_data$meta[train.idx, ],
  geneExp = train_data$geneExp[, train.idx], 
  ab = train_data$ab[, train.idx], 
  cytokine = train_data$cytokine[, train.idx],
  cell_freq = train_data$cell_freq[, train.idx]
)
train_data2$ab[1:5, 1:13]
lapply(train_data2, dim)

# RL: This creates a test dataset (test_data) by subsetting the original train_data list for the specimens that belong to the test set (test.idx). test_data: A list structured similarly to train_data, but containing only the test samples.
test_data = list(
  meta = train_data$meta[test.idx, ],
  geneExp = train_data$geneExp[, test.idx], 
  ab = train_data$ab[, test.idx], 
  cytokine = train_data$cytokine[, test.idx],
  cell_freq = train_data$cell_freq[, test.idx]
)
test_data$ab[1:5, 1:13]
lapply(test_data, dim)

train_data = train_data2
#save(train_data, file='processed_datasets/training_datasets/train_2020+21_TRAIN.RData')

train_data = test_data

#save(train_data, file='processed_datasets/training_datasets/train_2020+21_TEST.RData')
