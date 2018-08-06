setwd("~/Git-Projects/Git-Research-Projects/gene-set-feature-scoring/")
source("helperFunctions.R")
source("genomicFeatureAssignment.R")

gene_list <- read.table("./resources/genes_hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_list <- gene_list[, c(7,8,9,3)]
colnames(gene_list) <- c("chrom", "start", "end", "name")

reference <- "hN30"
res_dir <- "output/FACETS_Reference_hN30_8_2_18_1/" # Determine FACETS reference to use _ hN30

target_samples <- load_samples(classes = c("N","T", "F", "M"), sampleList = "./resources/sampleList.csv")
target_samples <- target_samples[target_samples != reference]
target_samples <- target_samples[target_samples != "hN32"]
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")

training_set <- retrieveTrainingSet(target_samples, gene_list, impute = TRUE, reference, res_dir)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)


aucData <- readRDS("./resources/listSampleTESAUC.RDS")
matrix_training_set <- attachLabelsToSet(matrix_training_set = matrix_training_set, labelData = aucData)

setwd("~/Git-Projects/Git-Research-Projects/gene-set-feature-scoring/")
write.csv(matrix_training_set, file ="mlOutput/coreTrainingSet_8_5_2018_1.csv")
